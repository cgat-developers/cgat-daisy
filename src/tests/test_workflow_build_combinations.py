import pytest
from daisy.workflow import build_combinations


def test_one_option():
    assert build_combinations(
        {"option1": ["value1", "value2"]}) == \
        [{'option1': 'value1'},
         {'option1': 'value2'}]


def test_two_options():
    assert build_combinations(
        {'option1': ["value1", "value2"],
         'option2': 'valueA'}) == \
         [{'option2': 'valueA', 'option1': 'value1'},
          {'option2': 'valueA', 'option1': 'value2'}]


def test_two_options():
    assert build_combinations(
        {'option1': ["value1", "value2"],
         'option2': ["valueA", "valueB"]}) == \
         [{'option1': 'value1', 'option2': 'valueA'},
          {'option1': 'value1', 'option2': 'valueB'},
          {'option1': 'value2', 'option2': 'valueA'},
          {'option1': 'value2', 'option2': 'valueB'}]


def test_complex_values():
    assert build_combinations(
        {'option1': [{"value1": [1, 2, 3]},
                     {"value2": [4, 5, 6]}]}) == \
        [{'option1': {'value1': [1, 2, 3]}},
         {'option1': {'value2': [4, 5, 6]}}]


def test_groupby_design(tmp_path):

    design_file = tmp_path / "design.tsv"
    with open(design_file, "w") as outf:
        outf.write("label\tc_option1\tc_option2\n")
        outf.write("label1\tvalue1\tvalueA\n")
        outf.write("label2\tvalue1\tvalueB\n")
        outf.write("label3\tvalue2\tvalueA\n")
        outf.write("label4\tvalue2\tvalueB\n")

    assert build_combinations(
        {"groupby": "file",
         "label": "label",
         "input": design_file,
         "option1": "c_option1",
         "option2": "c_option2"}) == \
         [{'option1': 'value1', 'option2': 'valueA', "name": "label1"},
          {'option1': 'value1', 'option2': 'valueB', "name": "label2"},
          {'option1': 'value2', 'option2': 'valueA', "name": "label3"},
          {'option1': 'value2', 'option2': 'valueB', "name": "label4"}]


def test_groupby_design_with_constant_option(tmp_path):

    design_file = tmp_path / "design.tsv"
    with open(design_file, "w") as outf:
        outf.write("label\tc_option1\tc_option2\n")
        outf.write("label1\tvalue1\tvalueA\n")
        outf.write("label2\tvalue1\tvalueB\n")
        outf.write("label3\tvalue2\tvalueA\n")
        outf.write("label4\tvalue2\tvalueB\n")

    assert build_combinations(
        {"groupby": "file",
         "label": "label",
         "input": design_file,
         "option1": "c_option1",
         "option2": "c_option2",
         "option3": "valueX"}) == \
         [{'option1': 'value1', 'option2': 'valueA', "name": "label1", "option3": "valueX"},
          {'option1': 'value1', 'option2': 'valueB', "name": "label2", "option3": "valueX"},
          {'option1': 'value2', 'option2': 'valueA', "name": "label3", "option3": "valueX"},
          {'option1': 'value2', 'option2': 'valueB', "name": "label4", "option3": "valueX"}]


def test_groupby_design_with_combinatorial_option(tmp_path):

    design_file = tmp_path / "design.tsv"
    with open(design_file, "w") as outf:
        outf.write("label\tc_option1\tc_option2\n")
        outf.write("label1\tvalue1\tvalueA\n")
        outf.write("label2\tvalue1\tvalueB\n")
        outf.write("label3\tvalue2\tvalueA\n")
        outf.write("label4\tvalue2\tvalueB\n")
    
    assert build_combinations(
        {"groupby": "file",
         "label": "label",
         "input": design_file,
         "option1": "c_option1",
         "option2": "c_option2",
         "option3": ["valueX", "valueY"]}) == \
         [{'option1': 'value1', 'option2': 'valueA', "name": "label1", "option3": "valueX"},
          {'option1': 'value1', 'option2': 'valueA', "name": "label1", "option3": "valueY"},
          {'option1': 'value1', 'option2': 'valueB', "name": "label2", "option3": "valueX"},
          {'option1': 'value1', 'option2': 'valueB', "name": "label2", "option3": "valueY"},
          {'option1': 'value2', 'option2': 'valueA', "name": "label3", "option3": "valueX"},
          {'option1': 'value2', 'option2': 'valueA', "name": "label3", "option3": "valueY"},
          {'option1': 'value2', 'option2': 'valueB', "name": "label4", "option3": "valueX"},
          {'option1': 'value2', 'option2': 'valueB', "name": "label4", "option3": "valueY"}]


def test_groupby_regex(tmp_path):

    assert build_combinations(
        {"groupby": "regex",
         "files_a": ["{}/data_0.a".format(tmp_path),
                     "{}/data_1.a".format(tmp_path)],
         "files_b": ["{}/data_0.b".format(tmp_path),
                     "{}/data_1.b".format(tmp_path)],
         "files_a_regex": r"data_(\d+).a",
         "files_b_regex": r"data_(\d+).b"}) == \
        [{'files_a': "{}/data_0.a".format(tmp_path),
          'files_b': "{}/data_0.b".format(tmp_path),
          'name': "0"},
         {'files_a': "{}/data_1.a".format(tmp_path),
          'files_b': "{}/data_1.b".format(tmp_path),
          'name': "1"}]


def test_groupby_regex_filters_when_data_point_missing(tmp_path):
    assert build_combinations(
        {"groupby": "regex",
         "files_a": ["{}/data_0.a".format(tmp_path)],
         "files_b": ["{}/data_0.b".format(tmp_path),
                     "{}/data_1.b".format(tmp_path)],
         "files_a_regex": r"data_(\d+).a",
         "files_b_regex": r"data_(\d+).b"}) == \
        [{'files_a': "{}/data_0.a".format(tmp_path),
          'files_b': "{}/data_0.b".format(tmp_path),
          'name': "0"}]


def test_groupby_regex_with_constant(tmp_path):

    assert build_combinations(
        {"groupby": "regex",
         "files_x": "x.y",
         "files_a": ["{}/data_0.a".format(tmp_path),
                     "{}/data_1.a".format(tmp_path)],
         "files_b": ["{}/data_0.b".format(tmp_path),
                     "{}/data_1.b".format(tmp_path)],
         "files_a_regex": r"data_(\d+).a",
         "files_b_regex": r"data_(\d+).b"}) == \
        [
            {'files_a': "{}/data_0.a".format(tmp_path),
             'files_b': "{}/data_0.b".format(tmp_path),
             'files_x': "x.y",
             'name': "0"},
            {'files_a': "{}/data_1.a".format(tmp_path),
             'files_b': "{}/data_1.b".format(tmp_path),
             'files_x': "x.y",
             'name': "1"},
        ]


def test_groupby_regex_with_combinatorial_option(tmp_path):

    assert build_combinations(
        {"groupby": "regex",
         "files_x": ["y.x", "z.x"],
         "files_a": ["{}/data_0.a".format(tmp_path),
                     "{}/data_1.a".format(tmp_path)],
         "files_b": ["{}/data_0.b".format(tmp_path),
                     "{}/data_1.b".format(tmp_path)],
         "files_a_regex": r"data_(\d+).a",
         "files_b_regex": r"data_(\d+).b"}) == \
        [
            {'files_a': "{}/data_0.a".format(tmp_path),
             'files_b': "{}/data_0.b".format(tmp_path),
             'files_x': "y.x",
             'name': "0"},
            {'files_a': "{}/data_0.a".format(tmp_path),
             'files_b': "{}/data_0.b".format(tmp_path),
             'files_x': "z.x",
             'name': "0"},
            {'files_a': "{}/data_1.a".format(tmp_path),
             'files_b': "{}/data_1.b".format(tmp_path),
             'files_x': "y.x",
             'name': "1"},
            {'files_a': "{}/data_1.a".format(tmp_path),
             'files_b': "{}/data_1.b".format(tmp_path),
             'files_x': "z.x",
             'name': "1"},
        ]


def test_groupby_named_regex(tmp_path):

    assert build_combinations(
        {"groupby": "regex",
         "files_a": ["{}/data_0.a".format(tmp_path),
                     "{}/data_1.a".format(tmp_path)],
         "files_b": ["{}/data_0.b".format(tmp_path),
                     "{}/data_1.b".format(tmp_path)],
         "files_a_regex": r"data_(?P<key1>\d+).a",
         "files_b_regex": r"data_(?P<key1>\d+).b"}) == \
        [{'files_a': "{}/data_0.a".format(tmp_path),
          'files_b': "{}/data_0.b".format(tmp_path),
          'name': "0"},
         {'files_a': "{}/data_1.a".format(tmp_path),
          'files_b': "{}/data_1.b".format(tmp_path),
          'name': "1"}]


def test_groupby_named_regex_paired(tmp_path):

    assert build_combinations(
        {"groupby": "regex",
         "files_a": ["{}/data_0_2.a".format(tmp_path),
                     "{}/data_0_3.a".format(tmp_path),
                     "{}/data_1_2.a".format(tmp_path),
                     "{}/data_1_3.a".format(tmp_path)],
         "files_b": ["{}/data_0.b".format(tmp_path),
                     "{}/data_1.b".format(tmp_path)],
         "files_a_regex": r"data_(?P<key1>\d+)_(?P<key2>\d+).a",
         "files_b_regex": r"data_(?P<key1>\d+).b"}) == \
        [{'files_a': "{}/data_0_2.a".format(tmp_path),
          'files_b': "{}/data_0.b".format(tmp_path),
          'name': "0_2"},
         {'files_a': "{}/data_0_3.a".format(tmp_path),
          'files_b': "{}/data_0.b".format(tmp_path),
          'name': "0_3"},
         {'files_a': "{}/data_1_2.a".format(tmp_path),
          'files_b': "{}/data_1.b".format(tmp_path),
          'name': "1_2"},
         {'files_a': "{}/data_1_3.a".format(tmp_path),
          'files_b': "{}/data_1.b".format(tmp_path),
          'name': "1_3"}]
