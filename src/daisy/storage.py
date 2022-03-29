"""Storage.py - Uploading metrics to storage
============================================

The :mod:`Storage` module contains various utility functions to upload
metrics computed by the daisy :ref:` `Benchmark` system into an SQL
database.

The main entry point to functionality in this module is the function
:meth:`upload_result`.

API
---

"""

import multiprocessing
from multiprocessing.util import Finalize
import os
import sys
import json
import re
import glob
import datetime
import pathlib
import pandas
import pandas.io.sql
import tqdm

import sqlalchemy
from sqlalchemy import inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import Column, Integer, String, DateTime, ForeignKey, text, \
    PrimaryKeyConstraint, LargeBinary

from sqlalchemy.dialects import postgresql
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import OperationalError

import cgatcore.pipeline as P
import cgatcore.iotools as IOTools
import cgatcore.experiment as E

from daisy.toolkit import touch, read_data, hash
from daisy.table_cache import TableCache, create_engine


######################################################
# Database schema

# Think about how best to use postgres schemas
# Currently there is one 'run' table across all
# to keep ids consistent
Base = declarative_base()

IS_POSTGRES = False
if IS_POSTGRES:
    JSONType = postgresql.JSON
else:
    JSONType = String


class BenchmarkRun(Base):
    __tablename__ = 'run'

    id = Column(Integer, primary_key=True)
    author = Column(String)
    created = Column(DateTime)
    pipeline_name = Column(String)
    pipeline_version = Column(String)
    pipeline_dir = Column(String)
    config = Column(JSONType)
    config_hash = Column(String)
    title = Column(String)
    description = Column(String)
    status = Column(String)

    # AH: for full ORM, activate references, but
    # not sure how this plays with the metric tables who
    # are created outside the ORM.
    # instances = relationship("BenchmarkInstance",
    #                          order_by="BenchmarkInstance.id",
    #                          backref="run")


class BenchmarkInstance(Base):
    __tablename__ = 'instance'

    id = Column(Integer, primary_key=True)
    run_id = Column(Integer, ForeignKey('run.id'))
    replication_id = Column(Integer)
    completed = Column(DateTime)
    input = Column(String)
    input_alias = Column(String)
    metric_name = Column(String)
    metric_version = Column(String)
    metric_options = Column(String)
    metric_hash = Column(String)
    metric_alias = Column(String)
    tool_name = Column(String)
    tool_version = Column(String)
    tool_options = Column(String)
    tool_hash = Column(String)
    tool_alias = Column(String)
    split_name = Column(String)
    split_version = Column(String)
    split_options = Column(String)
    split_hash = Column(String)
    split_alias = Column(String)
    split_subset = Column(String)
    meta_data = Column(JSONType)

    # run = relationship("BenchmarkRun",
    #                    backref='instances',
    #                    order_by=id)


# todo: add index on tag
# todo: separate table for tag and tags (normalized design)
class BenchmarkTag(Base):
    __tablename__ = 'tags'

    run_id = Column(Integer, ForeignKey('run.id'))
    tag = Column(String)

    __table_args__ = (PrimaryKeyConstraint(run_id, tag),)


class BenchmarkMetricStorage(Base):
    __tablename__ = 'metric_storage'

    run_id = Column(Integer, ForeignKey('run.id'))
    tablename = Column(String, index=True)
    schema = Column(String, index=True)
    bytes_uploaded = Column(Integer)
    bytes_resident = Column(Integer)

    __table_args__ = (PrimaryKeyConstraint(run_id, tablename),)


class BenchmarkBinaryData(Base):
    __tablename__ = 'binary_data'

    instance_id = Column(Integer, ForeignKey('instance.id'))
    filename = Column(String)
    path = Column(String)
    data = Column(LargeBinary)

    __table_args__ = (PrimaryKeyConstraint(instance_id, filename),)


def create_database(engine):

    # create_all only creates table that are not present
    # so it is save to call repeatedly.
    Base.metadata.create_all(engine)
    session = sessionmaker(bind=engine)()
    session.commit()


def save_benchmark_timings(path, tablename, table_cache, instance_id: int):

    fn = os.path.join(path, "benchmark.bench")
    if not os.path.exists(fn):
        P.get_logger().warn(
            "file {} does not exist, no tool timings uploaded".format(
                fn))
    else:
        table = pandas.read_csv(fn, sep="\t")
        table["instance_id"] = instance_id
        table_cache.add_table(table, tablename)


def transform_table_before_upload(tablename, table, instance_id: int, meta_data, table_cache):

    dtypes = None

    logger = P.get_logger()

    # melt table if set by metric
    if "metric_upload_melted" in meta_data:
        if tablename in meta_data["metric_upload_melted"]:
            melt_data = meta_data["metric_upload_melted"][tablename]
            table = pandas.melt(
                table,
                id_vars=melt_data.get("id_vars", None),
                value_vars=melt_data.get("value_vars", None),
                var_name=melt_data.get("var_name", None),
                value_name=melt_data.get("value_name", None))
            logger.debug("melted data from table {}".format(tablename))

    if "metric_upload_transpose" in meta_data:
        if tablename in meta_data["metric_upload_transpose"]:
            table = table.transpose()

    # upload into a separate table suffixed by instance id
    if "metric_upload_separate" in meta_data:
        if tablename in meta_data["metric_upload_separate"]:
            tablename = "{}_{}".format(tablename, instance_id)

    # normalize table by factorizing a column and storing its ids
    # in a separate table
    if "metric_upload_normalize" in meta_data:
        if tablename in meta_data["metric_upload_normalize"]:
            for column in meta_data["metric_upload_normalize"][tablename]:
                if column not in table.columns:
                    raise ValueError(
                        "unknown column {} in table {} to be normalized".format(
                            tablename, column))
                factors, names = table[column].factorize()
                table[column] = factors
                table.rename(columns={column: column + "_id"}, inplace=True)

                factor_table = pandas.DataFrame(
                    {column: names,
                     "id": list(range(len(names)))})
                factor_table["instance_id"] = instance_id
                table_cache.add_table(factor_table, tablename + "_factors")

    # store table as a matrix
    if "metric_upload_as_matrix" in meta_data:
        if tablename in meta_data["metric_upload_as_matrix"]:
            groupby_columns = meta_data["metric_upload_as_matrix"][tablename]
            if not isinstance(groupby_columns, list):
                groupby_columns = [groupby_columns]
            take_columns = [x for x in table.columns if x not in groupby_columns]
            row_index_column = take_columns.pop(0)
            rows = []
            if not groupby_columns:
                matrix = table.as_matrix(take_columns)
                rows.append([",".join(map(str, table[row_index_column])),
                             ",".join(map(str, take_columns)),
                             str(matrix.dtype),
                             matrix.tostring()])
            else:
                for key, group in table.groupby(by=groupby_columns):
                    if not isinstance(key, tuple):
                        key = [key]
                    matrix = group.as_matrix(take_columns)
                    rows.append(list(key) +
                                [",".join(map(str, group[row_index_column])),
                                 ",".join(map(str, take_columns)),
                                 str(matrix.dtype),
                                 matrix.tostring()])
            table = pandas.DataFrame.from_records(
                rows,
                columns=groupby_columns + ["rows", "columns", "dtype", "data"])
            dtypes = {"data": LargeBinary}

    return tablename, table, dtypes


def divine_paths(infile):
    """find tool/metric paths for infile"""

    # walk up the path to find "benchmark.info" as it might be
    # located on a higher level if the tool output multiple files.
    parts = list(pathlib.Path(os.path.dirname(infile)).parts)
    info_paths = []
    rootdir = os.getcwd()
    while parts:
        p = os.path.join(*parts)
        if p == rootdir:
            break
        if os.path.exists(os.path.join(p, "benchmark.info")):
            info_paths.append(p)
        parts.pop()
    info_paths = info_paths[::-1]

    # the level of nesting determines the layout:
    # 1 level: aggregation: tool == metric
    # 2 levels: tool + metric
    # 3 levels: tool + split + metric
    if len(info_paths) not in (1, 2, 3):
        raise ValueError(
            "for {}, expected two or three paths with info, "
            "got {}".format(infile, len(info_paths)))

    meta_data = {}

    if len(info_paths) == 1:
        tool_dir = metric_dir = info_paths[0]
        split_dir = None
    elif len(info_paths) == 2:
        tool_dir, metric_dir = info_paths
        split_dir = None
        # If there are multiple output files in aggregation, use
        # intermediate paths as split_subset factors.
        td = len(tool_dir.split(os.sep))
        tm = len(metric_dir.split(os.sep))
        delta = tm - td
        if delta > 1:
            meta_data["split_subset"] = re.sub(
                ".dir", "", os.sep.join(metric_dir.split(os.sep)[td:-1]))
    elif len(info_paths) == 3:
        tool_dir, split_dir, metric_dir = info_paths

    if tool_dir:
        md = read_data(os.path.join(tool_dir, "benchmark.info"), prefix="tool_")
        if "tool_action" in md and md["tool_action"] != "tool":
            raise ValueError("action for tool info {} is not 'tool', but '{}'".format(
                os.path.join(metric_dir, "benchmark.info"),
                md["tool_action"]))
        meta_data.update(md)

    if metric_dir:
        md = read_data(os.path.join(metric_dir, "benchmark.info"),
                       prefix="metric_")
        if "metric_action" in md:
            # ignore splits, they will be added through metrics
            if md["metric_action"] == "split":
                return None, None, None
            if md["metric_action"] != "metric":
                return tool_dir, None, None
        meta_data.update(md)

    if split_dir:
        md = read_data(os.path.join(split_dir, "benchmark.info"),
                       prefix="split_")
        if "split_action" in md and md["split_action"] != "split":
            raise ValueError("action for split info {} is not 'split', but '{}'".format(
                os.path.join(metric_dir, "benchmark.info"),
                md["split_action"]))

        meta_data.update(md)
        subset = os.path.basename(os.path.dirname(info_paths[-1]))
        if subset.endswith(".dir"):
            subset = subset[:-len(".dir")]
        meta_data["split_subset"] = subset

    return tool_dir, metric_dir, meta_data


def save_metric_data(meta_data, table_cache, schema, instance_id: int, session):

    logger = P.get_logger()
    metric_table_filter = None
    if "metric_no_upload" in meta_data:
        if meta_data["metric_no_upload"] == "*":
            logger.warn("upload turned off for metric {}".format(
                meta_data["metric_name"]))
            return
        else:
            metric_table_filter = re.compile(meta_data["metric_no_upload"])

    # multiple tablenames for multiple metric output
    #
    # Tables are added into schemas to avoid cluttering
    # the public namespace.
    # (if only blobs, no metric output file)
    if "metric_output_files" in meta_data:
        assert len(meta_data["metric_output_files"]) == \
            len(meta_data["metric_tablenames"])

        for output_file, tablename in zip(
                meta_data["metric_output_files"],
                meta_data["metric_tablenames"]):

            if metric_table_filter and metric_table_filter.search(tablename):
                logger.warn("upload for table {} turned off".format(
                    tablename))
                continue

            if not os.path.exists(output_file):
                logger.warning("output file {} does not exist - ignored".format(
                    output_file))
                continue

            if IOTools.is_empty(output_file):
                logger.warn("output file {} is empty - ignored".format(
                    output_file))
                continue

            # table = pandas.DataFrame({"values": [1, 2]})
            try:
                table = pandas.read_csv(output_file,
                                        sep="\t",
                                        comment="#",
                                        skip_blank_lines=True)
            except ValueError as e:
                logger.warn("table {} can not be read: {}".format(
                    output_file, str(e)))
                continue
            except pandas.parser.CParserError as e:
                logger.warn("malformatted table {} can not be read: {}".format(
                    output_file, str(e)))
                continue

            if table.empty:
                logger.warn("table {} is empty - ignored".format(output_file))
                continue

            tablename, table, dtypes = transform_table_before_upload(tablename,
                                                                     table,
                                                                     instance_id,
                                                                     meta_data,
                                                                     table_cache)

            if schema is None:
                tn = tablename
            else:
                tn = "{}.{}".format(schema, tablename)

            # add foreign key
            table["instance_id"] = instance_id
            logger.debug(f"saving data {table.shape} from {output_file} to table {tn} under {instance_id}")
            table_cache.add_table(table, tablename, dtypes)

    if "metric_blob_globs" in meta_data:
        metric_dir = meta_data["metric_outdir"]
        files = [glob.glob(os.path.join(metric_dir, x))
                 for x in meta_data["metric_blob_globs"]]
        files = IOTools.flatten(files)
        logger.debug(
            "uploading binary data in {} files from {} to "
            "table binary_data".format(len(files), metric_dir))
        table = []
        for fn in files:
            with IOTools.open_file(fn, "rb", encoding=None) as inf:
                data_row = BenchmarkBinaryData(
                    instance_id=instance_id,
                    filename=os.path.basename(fn),
                    path=fn,
                    data=inf.read())
                session.add(data_row)
            session.commit()

    if meta_data.get("metric_tableindices", None):
        table_cache.add_indices(meta_data["metric_tableindices"])


def instantiate_metrics(metrics, session, run_id):

    tool_dirs = set()
    for tool_dir, metric_dir, meta_data, mtime in metrics:

        if meta_data is None:
            continue

        # tool_input_files can either be a dictionary if a tool
        # or a simple list if aggregation.
        try:
            tool_input_files = [x["path"] for x in meta_data["tool_input_files"]]
        except TypeError:
            tool_input_files = meta_data["tool_input_files"]

        try:
            instance = BenchmarkInstance(
                run_id=run_id,
                replication_id=meta_data.get("tool_replication_id", 1),
                completed=datetime.datetime.fromtimestamp(mtime),
                input=",".join(tool_input_files),
                input_alias=meta_data["tool_input_alias"],
                tool_name=meta_data["tool_name"],
                tool_version=meta_data["tool_version"],
                tool_options=meta_data["tool_options"],
                tool_hash=meta_data["tool_option_hash"],
                tool_alias=meta_data.get("tool_alias", ""),
                metric_name=meta_data["metric_name"],
                metric_version=meta_data["metric_version"],
                metric_options=meta_data["metric_options"],
                metric_hash=meta_data["metric_option_hash"],
                metric_alias=meta_data.get("metric_alias", ""),
                split_name=meta_data.get("split_name", ""),
                split_version=meta_data.get("split_version", ""),
                split_options=meta_data.get("split_options", ""),
                split_hash=meta_data.get("split_option_hash", ""),
                split_alias=meta_data.get("split_alias", ""),
                split_subset=meta_data.get("split_subset", "all"),
                meta_data=json.dumps(meta_data))
        except KeyError as e:
            raise KeyError("missing required attribute {} in meta_data: {}".format(
                str(e), str(meta_data)))

        session.add(instance)
        session.commit()

        upload_tool_metrics = True
        if tool_dir:
            upload_tool_metrics = tool_dir not in tool_dirs
            if upload_tool_metrics:
                tool_dirs.add(tool_dir)

        assert instance.id is not None
        yield(tool_dir, metric_dir, meta_data, instance.id, upload_tool_metrics)


def generate_metric(infile):

    tool_dir, metric_dir, meta_data = divine_paths(infile)

    return tool_dir, metric_dir, meta_data, os.path.getmtime(infile)


# global variables for multiprocessing
resource = None


class Resource(object):
    def __init__(self, database_url: str, schema: str):
        self.database_url = database_url
        self.schema = schema

    def __enter__(self):
        table_cache = TableCache(self.database_url, self.schema)
        E.debug(f"{os.getpid()}: created resource={id(self)}: cache={id(table_cache)}")
        self.table_cache = table_cache
        return self

    def __exit__(self, *args, **kwargs):
        E.debug(f"{os.getpid()}: resource={id(self)}: final table cache flush for cache={id(self.table_cache)} started")
        E.debug(f"{os.getpid()}: resource={id(self)}: cache={id(self.table_cache)}: "
                f"{self.table_cache.get_cache_stats()}")
        self.table_cache.flush_all()
        E.debug(f"{os.getpid()}: resource={id(self)}: final table cache flush "
                f"for cache={id(self.table_cache)} completed")


def open_resource(url, schema):
    return Resource(url, schema)


def upload_metric(args):
    tool_dir, metric_dir, meta_data, instance_id, upload_tool_metrics = args

    if upload_tool_metrics:
        save_benchmark_timings(tool_dir,
                               "tool_timings",
                               upload_metric.table_cache,
                               instance_id)

    save_benchmark_timings(metric_dir,
                           "metric_timings",
                           upload_metric.table_cache,
                           instance_id)

    E.debug(f"{os.getpid()}: adding metrics for instance_id {instance_id} to cache={id(upload_metric.table_cache)}")
    save_metric_data(meta_data,
                     upload_metric.table_cache,
                     upload_metric.schema,
                     instance_id,
                     upload_metric.session)


def setup_worker(f, *args):

    global resource

    url, schema = args

    f.schema = schema
    Session = sessionmaker()
    f.session = Session()

    resource_cm = open_resource(url, schema)
    E.debug(f"{os.getpid()}: setting up worker for resource={id(resource)}")
    old_resource = resource
    resource = resource_cm.__enter__()
    E.debug(f"{os.getpid()}: new worker for resource={id(resource)} (old_resource={id(old_resource)})")

    # Register a finalizer to flush table cache
    Finalize(resource, resource.__exit__, exitpriority=16)

    E.debug(f"{os.getpid()}: adding cache={id(resource.table_cache)} from resource={id(resource)} "
            f"to worker={id(f)}, session={id(f.session)}")
    f.table_cache = resource.table_cache


def upload_metrics_tables(infiles: list,
                          run_id: int,
                          schema,
                          url: str,
                          max_workers: int = 10):

    logger = P.get_logger()

    engine = create_engine(url)

    Session = sessionmaker(bind=engine)
    session = Session()
    
    logger.info(f"{os.getpid()}: collecting upload items for {len(infiles)} input files")
    metric_f = generate_metric
    pool = multiprocessing.Pool(max_workers)
    metrics = pool.map(metric_f, infiles)
    pool.close()
    pool.join()

    logger.info(f"{os.getpid()}: instantiating {len(metrics)} metrics")
    data = list(tqdm.tqdm(instantiate_metrics(metrics, session, run_id),
                          total=len(metrics)))

    logger.info(f"{os.getpid()}: uploading {len(data)} items")
    upload_f = upload_metric
    initargs = (upload_f, url, schema)
    if max_workers == 1:
        setup_worker(*initargs)
        result = list(map(upload_f, data))
        global resource
        resource.table_cache.flush_all()
    else:
        logger.info(f"{os.getpid()}: loading data with {max_workers} cores")
        pool = multiprocessing.Pool(max_workers, initializer=setup_worker, initargs=initargs)
        pool.map(upload_f, data)
        pool.close()
        pool.join()


def mark_upload_complete(url: str, run_id: int):
    engine = create_engine(url)
    session = sessionmaker(bind=engine)()
    benchmark_run = session.query(BenchmarkRun).filter_by(id=run_id).first()
    benchmark_run.status = "complete"
    session.commit()


def upload_result(infiles, outfile, *extras):
    """upload results into database.

    Connection details for the database are taken from the
    configuration dictionary given as first argument to extras.  The
    configuration directory should have an element 'database' with the
    required field ``url`` and the optional field ``schema``.  For
    example, to upload to an sqlite database in the current directory
    called csvdb, use::

        config = {"database": {"url": "sqlite:///./csvdb"}}

    To use multiple cores, try::

        config = {"database": {"url": "sqlite:///./csvdb", "cores": 10}}

    Arguments
    ---------
    infiles: list
       List of files to upload. These should be the output
       of metric tasks in a benchmarking workflow.
    outfile: output file
       On success, an empty output file is created.
    extras: list
       List of one element containing a configuration directory
       (see above).

    """

    logger = P.get_logger()

    if len(extras) != 1:
        raise ValueError(
            "expecting only one extra argument "
            "(configuration dictionary)")

    config = extras[0]

    url = config["database"]["url"]
    max_workers = config["database"].get("cores", 1)
    
    schema = config["database"].get("schema", None)
    # TODO: check if schema exists to avoid incomplete
    # transaction.

    engine = create_engine(url)

    # Catch exceptions until database access on thame available
    try:
        create_database(engine)
    except OperationalError as msg:
        logger.warn(
            "could not connect to database at {}. "
            "The data will not be uploaded. Msg={}".format(
                url, str(msg)))
        return

    # Create schema if not exists
    if schema is not None:
        engine.execute(
            str(text("CREATE SCHEMA IF NOT EXISTS {}".format(schema))))

    pipeline_name = os.path.basename(sys.argv[0])
    logger.debug("uploading data to {}, schema={}".format(url, schema))
    # TODO: add dependencies
    # dependencies = infiles[1:]
    # meta_data = dict([("dependency{}".format(x), y) \
    #                  for x, y in enumerate(dependencies)])

    # need to set created dir somehow, important when re-loading
    # as otherwise all times will be the same.
    if os.path.exists("benchmark.yml"):
        s = os.stat("benchmark.yml")
        created = datetime.datetime.fromtimestamp(s.st_mtime)
    else:
        created = datetime.datetime.now()

    Session = sessionmaker(bind=engine)
    session = Session()

    benchmark_run = BenchmarkRun(
        author=os.environ.get("USER", "unknown"),
        # needs refactoring, should be: uploaded_at, created_at, run_at
        # uploaded_at=datetime.datetime.now(),
        created=created,
        pipeline_name=pipeline_name,
        pipeline_version=P.get_version().version,
        pipeline_dir=os.getcwd(),
        title=config["title"],
        description=config["description"],
        config=json.dumps(config),
        config_hash=hash(json.dumps(config)),
        status="incomplete")

    session.add(benchmark_run)
    session.commit()

    for tag in config["tags"]:
        benchmark_tag = BenchmarkTag(run_id=benchmark_run.id, tag=tag)
        session.add(benchmark_tag)

    session.commit()
    engine.dispose()
    del engine

    upload_metrics_tables(infiles,
                          benchmark_run.id,
                          schema,
                          url,
                          max_workers=max_workers)

    # upload table sizes
    # df_sizes = pandas.DataFrame.from_records(list(table_cache.uploaded_sizes.items()),
    #                                          columns=["tablename", "bytes_uploaded"])
    # df_sizes["bytes_resident"] = df_sizes.bytes_uploaded
    # df_sizes["run_id"] = benchmark_run.id
    # df_sizes["schema"] = schema
    # save_table(df_sizes,
    #            engine,
    #            "metric_storage",
    #            schema=None,
    #            is_sqlite3=is_sqlite3)

    mark_upload_complete(url, benchmark_run.id)

    logger.info("uploaded results under run_id {}".format(benchmark_run.id))
    touch(outfile)


def export_result(infile, outfile, *extras):

    def _link(infile, outfile):
        dirname = os.path.dirname(outfile)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        os.link(infile, outfile)

    suffixes = [".tbi", ".bai"]
    for suffix in suffixes:
        for fn in glob.glob(infile + suffix):
            _link(fn, outfile + suffix)

    _link(infile, outfile)


def get_instance_ids_for_run_id(run_id, engine):
    """return list of instance ids for a run id."""
    sql_statement = "SELECT id FROM instance WHERE run_id={}".format(run_id)
    sql_output = pandas.read_sql(sqlalchemy.text(sql_statement), engine)
    instance_id_list = set(sql_output['id'].tolist())
    return instance_id_list


def purge_run_id(run_id, url, dry_run=False, schemas=None):
    """remove a run from a database.
    """

    logger = P.get_logger()
    engine = sqlalchemy.create_engine(url)
    connection = engine.connect()

    # automap
    metadata = sqlalchemy.MetaData()
    metadata.reflect(engine)
    base = automap_base(metadata=metadata)
    base.prepare()

    if schemas is None:
        insp = inspect(engine)
        schemas = insp.get_schema_names()
        # note: default sqlite schema is "main"
        if 'public' in schemas:
            schemas.remove('public')
        if 'information_schema' in schemas:
            schemas.remove('information_schema')

    logger.debug("getting instance_id list of run_id={}".format(run_id))
    instance_ids = set(get_instance_ids_for_run_id(run_id, engine))
    logger.debug("found {} instances for run_id={}".format(len(instance_ids), run_id))
    non_metric_tables = ['run',
                         'instance',
                         'binary_data',
                         'metric_timings',
                         'tool_timings',
                         'metric_storage',
                         'tags']

    # delete from tables with field "instance_id"
    if instance_ids:
        for schema in schemas:
            # automap the schema
            metadata_schema = sqlalchemy.MetaData()
            metadata_schema.reflect(engine, schema=schema)
            base_schema = automap_base(metadata=metadata_schema)
            base_schema.prepare()
            for table_name in list(base_schema.metadata.tables.keys()):
                table = sqlalchemy.Table(table_name,
                                         metadata_schema,
                                         autoload=True)
                if "instance_id" not in table.c:
                    continue
                logger.debug("deleting data in {}".format(table_name))
                delete = table.delete().where(
                    table.c.instance_id.in_(instance_ids))
                if not dry_run:
                    connection.execute(delete)

    # delete from tables with field "run_id"
    for table_name in base.metadata.tables.keys():
        table = sqlalchemy.Table(table_name, metadata, autoload=True)
        if "run_id" not in table.c:
            continue
        logger.info("deleting data in {} for run_id {}".format(table_name, run_id))
        delete = table.delete().where(table.c.run_id == run_id)
        if not dry_run:
            connection.execute(delete)

    table = sqlalchemy.Table('run', metadata, autoload=True)
    delete = table.delete().where(table.c.id == run_id)
    logger.info("deleting data in 'run' for id {}".format(run_id))
    if not dry_run:
        connection.execute(delete)
