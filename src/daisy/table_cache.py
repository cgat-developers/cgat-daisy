import os
import re
import collections
import sqlite3
import time
import pandas
import sqlalchemy
from sqlalchemy import inspect
from sqlalchemy import text
import cgatcore.pipeline as P
import cgatcore.experiment as E


RESERVED_WORDS = {
    "all": "total"}


@E.cached_function
def get_columns(tablename, engine):
    """return list of column names in table"""
    metadata = sqlalchemy.MetaData(engine)
    tb = sqlalchemy.Table(tablename, metadata, autoload=True,
                          autoload_with=engine)
    return [x.name for x in tb.columns]


def create_engine(url: str):
    is_sqlite3 = url.startswith("sqlite")
    if is_sqlite3:
        connect_args = {'check_same_thread': False}
    else:
        connect_args = {}
    engine = sqlalchemy.create_engine(
        url,
        connect_args=connect_args)
    return engine


def sql_sanitize_columns(columns):

    # special chars
    columns = [re.sub(r"[\[\]().,:]", "_", str(x)) for x in columns]

    columns = [re.sub("%", "percent", x)
               for x in columns]

    columns = [x.lower() for x in columns]

    columns = [RESERVED_WORDS.get(x, x) for x in columns]
    return columns


class TableExistsException(Exception):
    pass


class TableSchemaChangedException(Exception):
    pass


class IndexExistsException(Exception):
    pass


def retry_table_to_sql(table, *args, **kwargs):
    """retry SQL statements retrying when database is locked"""
    while True:
        try:
            table.to_sql(*args, **kwargs)
        except (sqlalchemy.exc.OperationalError, sqlite3.OperationalError) as ex:
            msg = str(ex)
            if "database is locked" in msg:
                E.debug("database is locked")
                time.sleep(1)
                continue
            elif "database schema has changed" in msg:
                E.debug("schema has changed")
                time.sleep(1)
                continue
            elif "already exists" in msg:
                raise TableExistsException(str(ex))
            raise
        except pandas.io.sql.DatabaseError as ex:
            msg = str(ex)
            if "database schema has changed" in msg:
                E.debug("schema has changed")
                time.sleep(1)
                continue
                # raise TableSchemeChangedException(
                #     f"database schema changed: args={args}, kwargs={kwargs}")
            raise
        except ValueError as ex:
            # pandas throws ValueError
            if "already exists" in str(ex):
                raise TableExistsException(str(ex))
            raise
        except Exception as ex:
            raise
        break


def retry_sql_execute(engine, statement):
    while True:
        try:
            engine.execute(statement)
        except (sqlalchemy.exc.OperationalError, sqlite3.OperationalError) as ex:
            msg = str(ex)
            if "database is locked" in msg:
                E.debug("database is locked")
                time.sleep(1)
                continue
            elif "database schema has changed" in msg:
                E.debug("schema has changed")
                time.sleep(1)
                continue
            elif "index" in msg and "already exists" in msg:
                raise IndexExistsException(msg)
            else:
                raise
        except Exception as ex:
            raise


def add_columns_to_table(columns, table, tablename, engine):

    pandas_engine = pandas.io.sql.SQLDatabase(engine)
    pandas_table = pandas.io.sql.SQLTable(tablename,
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
            try:
                retry_sql_execute(engine, statement)
            except (sqlalchemy.exc.OperationalError, sqlite3.OperationalError) as ex:
                if "duplicate column name" not in str(ex):
                    raise


def reconcile_columns(tablename, engine, table):

    logger = P.get_logger()
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


def save_table(table: pandas.DataFrame,
               url: str,
               tablename: str,
               schema: str = None,
               dtypes=None,
               indices=["instance_id"]):
    logger = P.get_logger()
    table.columns = sql_sanitize_columns(table.columns)

    engine = create_engine(url)

    # pandas/sqlite3 prefers the raw connection, otherwise error:
    # AttributeError: 'Engine' object has no attribute 'rollback'
    if url.startswith("sqlite"):
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

    try:
        retry_table_to_sql(table,
                           tablename,
                           _engine,
                           schema=schema,
                           if_exists="fail",
                           index=False,
                           dtype=dtypes,
                           chunksize=sql_chunk_size)
        E.debug(f"table {tablename} was new")
        create_index = True
    except TableExistsException:
        E.debug(f"table {tablename} already exists - appending")

    if create_index:
        # sqlite requires an index name
        if schema:
            tablename = "{}.{}".format(schema, tablename)

        for field in indices:
            E.debug(f"creating index on {field} for {tablename}")
            try:
                retry_sql_execute(
                    _engine,
                    str(text("CREATE INDEX {} ON {} ({})".format(
                        re.sub("[-.]", "_", tablename) + "_" + field,
                        tablename,
                        field))))
            except IndexExistsException:
                pass
            except TypeError as ex:
                logger.warn("could not create index: {}".format(str(ex)))
            except sqlalchemy.exc.ProgrammingError as ex:
                logger.warn("could not create index: {}".format(str(ex)))
    else:
        reconcile_columns(tablename, engine, table)
        retry_table_to_sql(table,
                           tablename,
                           _engine,
                           schema=schema,
                           if_exists="append",
                           index=False,
                           dtype=dtypes,
                           chunksize=sql_chunk_size)


class TableCache():

    # 1 Gb
    max_total_bytes = 1e9
    # 50 Mb
    max_table_bytes = 5e7

    def __init__(self, database_url, schema):
        self.database_url = database_url
        self.schema = schema
        self.cache = {}

        self.total_size = 0
        self.sizes = collections.defaultdict(int)
        self.uploaded_sizes = collections.defaultdict(int)
        self.dtypes = {}
        self.logger = P.get_logger()
        self.indices = {}
        self.have_created_indices = False

    def flush_table(self, tablename):

        table = self.cache[tablename]

        self.logger.debug(f"{os.getpid()}: uploading table {tablename} from cache={id(self)} "
                          f"with {self.sizes[tablename]} bytes")
        save_table(table,
                   self.database_url,
                   tablename,
                   self.schema,
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

    def get_cache_stats(self):
        stats = {}
        for tablename in list(self.cache.keys()):
            stats[tablename] = {
                "rows": self.cache[tablename].shape[0],
                "columns": self.cache[tablename].shape[1]
            }
        return stats

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
            new_columns = set(table.columns) - set(self.cache[tablename].columns)
            if new_columns:
                self.logger.warning(f"additional columns for {tablename}: {new_columns}")
            try:
                self.cache[tablename] = pandas.concat([self.cache[tablename], table])
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

        if self.database_url is None:
            return

        try:
            engine = create_engine(self.database_url)
        except AttributeError:
            # do not create indices if sqlalchemy already shut down
            return

        if not self.have_created_indices:
            for index_name, info in self.indices.items():
                table_name, fields = info
                self.logger.debug("creating index {} on {} and fields {}".format(
                    index_name, table_name, fields))
                try:
                    engine.execute("CREATE INDEX {} ON {} ({})".format(
                        index_name, table_name, fields))
                except sqlalchemy.exc.OperationalError as ex:
                    self.logger.warn("could not create index: {}".format(str(ex)))
            self.have_created_indices = True

    def __del__(self):
        E.debug(f"closing table cache {id(self)}")
        self.close()
