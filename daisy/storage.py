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

import os
import sys
import json
import datetime
import collections
import pandas
import pandas.io.sql
import re
import glob

import sqlalchemy
from sqlalchemy.engine import reflection
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import Column, Integer, String, DateTime, ForeignKey, text, \
    PrimaryKeyConstraint, LargeBinary

from sqlalchemy.dialects import postgresql
from sqlalchemy.orm import sessionmaker
from sqlalchemy.exc import OperationalError

import cgatcore.pipeline as P
import cgatcore.iotools as IOTools

from daisy.toolkit import touch, read_data, hash
import cgatcore.experiment as E

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


RESERVED_WORDS = {
    "all": "total"}


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

    Session = sessionmaker(bind=engine)
    session = Session()
    session.commit()


@E.cached_function
def get_columns(tablename, engine):
    """return list of column names in table"""
    metadata = sqlalchemy.MetaData(engine)
    tb = sqlalchemy.Table(tablename, metadata, autoload=True,
                          autoload_with=engine)
    return [x.name for x in tb.columns]


def sql_sanitize_columns(columns):

    # special chars
    columns = [re.sub("[\[\]().,:]", "_", str(x)) for x in columns]

    columns = [re.sub("%", "percent", x)
               for x in columns]

    columns = [x.lower() for x in columns]

    columns = [RESERVED_WORDS.get(x, x) for x in columns]
    return columns


def add_columns_to_table(columns, table, tablename, engine):

    pandas_engine = pandas.io.sql.SQLDatabase(engine)
    pandas_table = pandas.io.sql.SQLTable(
        tablename,
        pandas_engine,
        frame=table)

    new_columns = set(columns)
    logger = P.get_logger()

    for column in pandas_table.table.columns:
        if column.name in new_columns:
            statement = "ALTER TABLE {} ADD COLUMN {} {}".format(
                tablename,
                column.name,
                column.type)

            logger.debug("SQL: {}".format(statement))
            engine.execute(statement)


def save_table(table, engine, tablename, schema=None,
               is_sqlite3=False,
               dtypes=None,
               indices=["instance_id"]):
    logger = P.get_logger()
    table.columns = sql_sanitize_columns(table.columns)

    # pandas/sqlite3 prefers the raw connection, otherwise error:
    # AttributeError: 'Engine' object has no attribute 'rollback'
    if is_sqlite3:
        _engine = engine.raw_connection()
        # In pandas >= 0.23 and using sqlite as a backend, the
        # pandas.DataFrame.to_sql command fails with "OperationalError:
        # (sqlite3.OperationalError) too many SQL variables". The reason is a
        # fixed limit in sqlite, SQLITE_MAX_VARIABLE_NUMBER, which is by
        # default set to 999.
        sql_chunk_size = 999 // (len(table.columns) + 1)
    else:
        _engine = engine
        sql_chunk_size = None

    # lower case all table names. Otherwise issues with psql
    # mixed case access
    tablename = tablename.lower()
    create_index = False

    if schema is not None:
        # for postgres tables within a schema, engine.has_table does
        # not work. A solution would be to set the search_path within
        # a connection. The code below fails if there is a column
        # mismatch.
        # TODO: investigate how to use a connection within pandas.
        try:
            # table is new, create index
            table.to_sql(tablename,
                         engine,
                         schema=schema,
                         if_exists="fail",
                         index=False,
                         dtype=dtypes,
                         chunksize=sql_chunk_size)
            create_index = True
        except ValueError:
            table.to_sql(tablename,
                         engine,
                         schema=schema,
                         if_exists="append",
                         index=False,
                         dtype=dtypes,
                         chunksize=sql_chunk_size)

    elif engine.has_table(tablename, schema):

        existing_columns = set(get_columns(tablename, engine))
        proposed_columns = set(table.columns)

        obsolete_columns = existing_columns.difference(proposed_columns)
        if obsolete_columns:
            logger.warn("the following columns are obsolete in {}: {}. "
                        "empty data will be inserted"
                        .format(tablename, ", ".join(obsolete_columns)))
            # create empty columns
            for column in obsolete_columns:
                table[column] = None

        new_columns = proposed_columns.difference(existing_columns)
        if new_columns:
            logger.warn("new columns found for {}: the following columns "
                        "will be added: {} ".format(
                            tablename,
                            ", ".join(new_columns)))

            add_columns_to_table(new_columns,
                                 table,
                                 tablename,
                                 engine)
            # clear cache of memoization function
            get_columns.delete(tablename, engine)

        # append
        table.to_sql(tablename,
                     engine,
                     schema=schema,
                     if_exists="append",
                     index=False,
                     dtype=dtypes,
                     chunksize=sql_chunk_size)
    else:
        # table is new, create index
        table.to_sql(tablename,
                     engine,
                     schema=schema,
                     if_exists="fail",
                     index=False,
                     dtype=dtypes,
                     chunksize=sql_chunk_size)
        create_index = True

    if create_index:
        # sqlite requires an index name
        if schema:
            tablename = "{}.{}".format(schema, tablename)

        for field in indices:
            try:
                engine.execute(
                    text("CREATE INDEX {} ON {} ({})".format(
                        re.sub("[-.]", "_", tablename) + "_" + field,
                        tablename,
                        field)))
            except sqlalchemy.exc.ProgrammingError as ex:
                logger.warn("could not create index: {}".format(str(ex)))


def save_benchmark_timings(path, tablename, engine, instance, schema,
                           is_sqlite3):

    fn = os.path.join(path, "benchmark.bench")
    if not os.path.exists(fn):
        P.get_logger().warn(
            "file {} does not exist, no tool timings uploaded".format(
                fn))
    else:
        try:
            tool_bench_data = pandas.read_csv(fn, sep="\t")
        except pandas.errors.EmptyDataError:
            P.get_logger().warn(f"file {fn} is empty, no tool timings uploaded")
            return
        tool_bench_data["instance_id"] = instance.id
        save_table(tool_bench_data, engine,
                   tablename,
                   schema=None,
                   is_sqlite3=is_sqlite3)


class TableCache():

    # 1 Gb
    max_total_bytes = 1e9
    # 50 Mb
    max_table_bytes = 5e7

    def __init__(self, engine, schema, is_sqlite3):
        self.engine = engine
        self.schema = schema
        self.cache = {}

        self.is_sqlite3 = is_sqlite3
        self.total_size = 0
        self.sizes = collections.defaultdict(int)
        self.uploaded_sizes = collections.defaultdict(int)
        self.dtypes = {}
        self.logger = P.get_logger()
        self.indices = {}
        self.have_created_indices = False

    def flush_table(self, tablename):

        table = self.cache[tablename]

        self.logger.debug("uploading table {}: {} bytes".format(
            tablename,
            self.sizes[tablename]))
        save_table(table,
                   self.engine,
                   tablename,
                   self.schema,
                   is_sqlite3=self.is_sqlite3,
                   dtypes=self.dtypes.get(tablename, None))
        del table
        del self.cache[tablename]
        self.uploaded_sizes[tablename] += self.sizes[tablename]
        self.sizes[tablename] = 0

    def flush_all(self):
        for tablename in list(self.cache.keys()):
            self.flush_table(tablename)

        self.total_size = 0
        self.cache = {}
        self.sizes = collections.defaultdict(int)

    def add_indices(self, indices):
        self.indices.update(indices)

    def add_table(self, table, tablename, dtypes=None):

        memory_usage = sum(table.memory_usage(deep=True))
        self.total_size += memory_usage
        self.sizes[tablename] += memory_usage
        self.dtypes[tablename] = dtypes

        if tablename not in self.cache:
            self.cache[tablename] = table
        else:
            if set(self.cache[tablename].columns) != set(table.columns):
                raise ValueError(
                    "column mismatch for {}: "
                    "table={}, cache={}".format(
                        tablename,
                        self.cache[tablename].columns,
                        table.columns))
            try:
                self.cache[tablename] = pandas.concat(
                    [self.cache[tablename],
                     table])
            except AssertionError:
                raise

        if self.total_size > self.max_total_bytes:
            self.logger.debug(
                "force full cache flush, memory usage = {}".format(
                    self.total_size))
            self.flush_all()

        if self.sizes[tablename] > self.max_table_bytes:
            self.logger.debug(
                "force cache flush for table {}, memory usage = {}".format(
                    tablename,
                    self.sizes[tablename]))

            self.flush_table(tablename)

    def close(self):
        self.flush_all()

        if not self.have_created_indices:
            for index_name, info in self.indices.items():
                table_name, fields = info
                self.logger.debug("creating index {} on {} and fields {}".format(
                    index_name, table_name, fields))
                try:
                    self.engine.execute("CREATE INDEX {} ON {} ({})".format(
                        index_name, table_name, fields))
                except sqlalchemy.exc.OperationalError as ex:
                    self.logger.warn("could not create index: {}".format(str(ex)))
            self.have_created_indices = True

    def __del__(self):
        self.close()


def transform_table_before_upload(tablename, table, instance, meta_data, table_cache):

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
            tablename = "{}_{}".format(tablename, instance.id)

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
                factor_table["instance_id"] = instance.id
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


def upload_result(infiles, outfile, *extras):
    """upload results into database.

    Connection details for the database are taken from the
    configuration dictionary given as first argument to extras.  The
    configuration directory should have an element 'database' with the
    required field ``url`` and the optional field ``schema``.  For
    example, to upload to an sqlite database in the current directory
    called csvdb, use::

        config = {"database": {"url": "sqlite:///./csvdb"}}

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
    is_sqlite3 = url.startswith("sqlite")

    if is_sqlite3:
        connect_args = {'check_same_thread': False}
    else:
        connect_args = {}

    schema = config["database"].get("schema", None)
    # TODO: check if schema exists to avoid incomplete
    # transaction.

    engine = sqlalchemy.create_engine(
        url,
        connect_args=connect_args)

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
            text("CREATE SCHEMA IF NOT EXISTS {}".format(schema)))

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

    Session = sessionmaker(bind=engine)
    session = Session()
    session.add(benchmark_run)
    session.commit()

    for tag in config["tags"]:
        benchmark_tag = BenchmarkTag(
            run_id=benchmark_run.id,
            tag=tag)
        session.add(benchmark_tag)
    session.commit()

    tool_dirs = set()

    table_cache = TableCache(engine, schema, is_sqlite3)

    for infile in infiles:

        path, name = os.path.split(infile)

        # walk up the path to find "benchmark.info" as it might be
        # located on a higher level if the tool output multiple files.
        parts = path.split(os.sep)

        info_paths = []
        rootdir = os.getcwd()
        while len(parts):
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
            d = tm - td
            if d > 1:
                meta_data["split_subset"] = re.sub(
                    ".dir", "",
                    os.sep.join(
                        metric_dir.split(os.sep)[td:-1]))
        elif len(info_paths) == 3:
            tool_dir, split_dir, metric_dir = info_paths

        if tool_dir:
            d = read_data(os.path.join(tool_dir, "benchmark.info"),
                          prefix="tool_")
            if "tool_action" in d:
                assert d["tool_action"] == "tool"
            meta_data.update(d)

        if metric_dir:
            d = read_data(os.path.join(metric_dir, "benchmark.info"),
                          prefix="metric_")
            if "metric_action" in d:
                # ignore splits, they will be added through metrics
                if d["metric_action"] == "split":
                    continue
                assert d["metric_action"] == "metric", \
                    "action for metric info {} is not 'metric', but '{}'" \
                    .format(os.path.join(metric_dir, "benchmark.info"),
                            d["metric_action"])

        meta_data.update(d)

        if split_dir:
            d = read_data(os.path.join(split_dir, "benchmark.info"),
                          prefix="split_")
            if "split_action" in d:
                assert d["split_action"] == "split"
            meta_data.update(d)
            subset = os.path.basename(
                os.path.dirname(info_paths[-1]))
            if subset.endswith(".dir"):
                subset = subset[:-len(".dir")]
            meta_data["split_subset"] = subset

        # tool_input_files can either be a dictionary if a tool
        # or a simple list if aggregation.
        try:
            tool_input_files = [x["path"] for x in meta_data["tool_input_files"]]
        except TypeError:
            tool_input_files = meta_data["tool_input_files"]

        try:
            instance = BenchmarkInstance(
                run_id=benchmark_run.id,
                replication_id=meta_data.get("tool_replication_id", 1),
                completed=datetime.datetime.fromtimestamp(
                    os.path.getmtime(infile)),
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
            raise KeyError("missing required attribute {} in {}".format(
                str(e), str(meta_data)))

        session.add(instance)
        session.commit()

        # avoid multiple upload of tool data
        if tool_dir and tool_dir not in tool_dirs:
            tool_dirs.add(tool_dir)
            save_benchmark_timings(tool_dir,
                                   "tool_timings",
                                   engine, instance, schema,
                                   is_sqlite3)

        save_benchmark_timings(metric_dir,
                               "metric_timings",
                               engine, instance, schema,
                               is_sqlite3)

        metric_table_filter = None
        if "metric_no_upload" in meta_data:
            if meta_data["metric_no_upload"] == "*":
                logger.warn("upload turned off for metric {}".format(
                    meta_data["metric_name"]))
                continue
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
                    logger.warn("output file {} does not exist - ignored".format(
                        output_file))
                    continue

                if IOTools.is_empty(output_file):
                    logger.warn("output file {} is empty - ignored".format(
                        output_file))
                    continue

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

                if len(table) == 0:
                    logger.warn("table {} is empty - ignored".format(output_file))
                    continue

                tablename, table, dtypes = transform_table_before_upload(tablename,
                                                                         table,
                                                                         instance,
                                                                         meta_data,
                                                                         table_cache)

                if schema is None:
                    tn = tablename
                else:
                    tn = "{}.{}".format(schema, tablename)

                logger.debug("saving data from {} to table {}".format(output_file, tn))
                # add foreign key
                table["instance_id"] = instance.id
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
                        instance_id=instance.id,
                        filename=os.path.basename(fn),
                        path=fn,
                        data=inf.read())
                    session.add(data_row)
                session.commit()

        if meta_data.get("metric_tableindices", None):
            table_cache.add_indices(meta_data["metric_tableindices"])

    table_cache.close()
    touch(outfile)

    # upload table sizes
    df_sizes = pandas.DataFrame.from_records(
        list(table_cache.uploaded_sizes.items()),
        columns=["tablename", "bytes_uploaded"])
    df_sizes["bytes_resident"] = df_sizes.bytes_uploaded
    df_sizes["run_id"] = benchmark_run.id
    df_sizes["schema"] = schema
    save_table(df_sizes,
               engine,
               "metric_storage",
               schema=None,
               is_sqlite3=is_sqlite3)

    benchmark_run.status = "complete"
    session.commit()

    engine.dispose()
    del engine

    logger.info("uploaded results under run_id {}".format(benchmark_run.id))


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
    engine = sqlalchemy.create_engine(url)
    connection = engine.connect()

    # automap
    metadata = sqlalchemy.MetaData()
    metadata.reflect(engine)
    base = automap_base(metadata=metadata)
    base.prepare()

    if schemas is None:
        insp = reflection.Inspector.from_engine(engine)
        schemas = insp.get_schema_names()
        # note: default sqlite schema is "main"
        if 'public' in schemas:
            schemas.remove('public')
        if 'information_schema' in schemas:
            schemas.remove('information_schema')

    E.debug("getting instance_id list of run_id={}".format(run_id))
    instance_ids = set(get_instance_ids_for_run_id(run_id, engine))
    E.debug("found {} instances for run_id={}".format(len(instance_ids), run_id))
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
                E.info("deleting data in {}".format(table_name))
                delete = table.delete().where(
                    table.c.instance_id.in_(instance_ids))
                # E.debug(delete)
                if not dry_run:
                    connection.execute(delete)

    # delete from tables with field "run_id"
    for table_name in base.metadata.tables.keys():
        table = sqlalchemy.Table(table_name, metadata, autoload=True)
        if "run_id" not in table.c:
            continue
        E.info("deleting data in {} for run_id {}".format(table_name, run_id))
        delete = table.delete().where(table.c.run_id == run_id)
        # E.debug(delete)
        if not dry_run:
            connection.execute(delete)

    table = sqlalchemy.Table('run', metadata, autoload=True)
    delete = table.delete().where(table.c.id == run_id)
    E.info("deleting data in 'run' for id {}".format(run_id))
    # E.debug(delete)
    if not dry_run:
        connection.execute(delete)
