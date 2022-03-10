import os
import json
import pandas
import pytest
import sqlalchemy
from daisy.storage import upload_result


@pytest.fixture
def benchmark_layout(tmp_path):

    ntools = 2
    nmetrics = 4000

    tools = [f"tool{x}" for x in range(ntools)]
    metrics = [f"metric{x}" for x in range(nmetrics)]

    def generate_bench(d):
        with open(os.path.join(d, "benchmark.bench"), "w") as outf:
            outf.write("statement\nmanual\n")

    def generate_info(d, name):
        with open(os.path.join(d, "benchmark.info"), "w") as outf:
            outf.write(json.dumps({
                "input_files": [],
                "version": "test_version",
                "options": "test_options",
                "option_hash": "test_option_hash",
                "name": name,
                "input_alias": "test_alias",
                "output_files": [os.path.join(d, f"{name}.tsv")],
                "tablenames": [f"{name}"],
            }))

    outfiles = []

    for tool in tools:
        tool_dir = os.path.join(tmp_path, f"{tool}.dir")
        os.makedirs(tool_dir)
        generate_bench(tool_dir)
        generate_info(tool_dir, name=tool)
        for metric in metrics:
            metric_dir = os.path.join(tool_dir, f"{metric}.dir")
            os.makedirs(metric_dir)
            generate_bench(metric_dir)
            generate_info(metric_dir, name=metric)
            metric_filename = os.path.join(metric_dir, f"{metric}.tsv")
            with open(metric_filename, "w") as outf:
                outf.write("metric_name\tmetric_value\n")
                outf.write(f"{metric}\t1\n")
            outfiles.append(metric_filename)
    return outfiles, tools, metrics


@pytest.mark.parametrize("max_workers", [1, 5])
def test_upload(benchmark_layout, tmp_path, max_workers):

    outfiles, tools, metrics = benchmark_layout
    db_path = f"sqlite:///{tmp_path}/csvdb"
    upload_result(outfiles,
                  os.path.join(tmp_path, "upload.log"),
                  {"title": "test",
                   "description": "test",
                   "tags": [],
                   "database": {"url": db_path}},
                  max_workers=max_workers)

    db_engine = sqlalchemy.create_engine(db_path)

    insp = sqlalchemy.inspect(db_engine)
    assert insp.get_table_names() == sorted([
        'binary_data', 'instance',
        'metric_storage', 'metric_timings', 'run',
        'tags', 'tool_timings'] + metrics)

    with db_engine.connect() as conn:
        run_df = pandas.read_sql("SELECT * FROM run", conn)
    assert len(run_df) == 1

    with db_engine.connect() as conn:
        ins_df = pandas.read_sql("SELECT * FROM instance", conn)
    assert list(ins_df.run_id.unique()) == [1]
    assert list(ins_df.id) == list(range(1, len(ins_df) + 1))

    for metric in metrics:
        with db_engine.connect() as conn:
            df = pandas.read_sql(f"SELECT * FROM {metric}", conn)
        assert len(df) == len(tools)

    with db_engine.connect() as conn:
        tt_df = pandas.read_sql("SELECT * FROM tool_timings", conn)
    assert len(tt_df) == len(tools)

    with db_engine.connect() as conn:
        mt_df = pandas.read_sql("SELECT * FROM metric_timings", conn)
    assert len(mt_df) == len(metrics) * len(tools)
