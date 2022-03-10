import re
import collections
import itertools
import pandas
import cgatcore.iotools


# parameters that will no be used to compute combinations
RESERVED_WORDS = ["ignore"]
PLACEHOLDER_KEY = "_placeholder"


def merge_shared_values(values):

    try:
        shared = [x for x in values if x.startswith("prefix")]
    except AttributeError:
        return values

    if shared:
        prefix = re.sub("prefix\s*=\s*", "", shared[0])
        return [prefix + " " + x for x in values if not x.startswith("prefix")]
    else:
        return values


def build_combinations_from_design_file(config):

    # use a design-file to define groups
    if "label" not in config:
        raise ValueError("using file requires a 'label' column to be set")

    label_columns = config["label"]
    if not isinstance(label_columns, list) or isinstance(label_columns, tuple):
        label_columns = [label_columns]

    filelist = config["input"]
    if not isinstance(filelist, list):
        filelist = [filelist]

    if len(filelist) > 1:
        raise NotImplementedError("using multiple files is not implemented")

    for fn in filelist:
        with cgatcore.iotools.open_file(fn) as inf:
            df = pandas.read_table(inf, dtype=str)

    for label_column in label_columns:
        if label_column not in df.columns:
            raise ValueError(
                "label column {} specified, but does not exist in {}".format(
                    label_column, fn))

    map_column2slot = {}
    shared_values = set()
    columns = set(df.columns)

    for key, value in list(config.items()):
        if key == "label":
            continue

        shared_value = True
        if not isinstance(value, list):
            value = [value]

        for v in value:
            if v in columns and v not in label_columns:
                map_column2slot[v] = key
                shared_value = False

        if shared_value and key != "input":
            shared_values.add(key)

    if len(map_column2slot) == 0:
        raise ValueError(
            "no mapping found between column headers ({}) "
            "and slots in config file ({})".format(
                ",".join(df.columns),
                ",".join(list(config.keys()))))

    combinations = []
    for row in df.iterrows():
        combination = {}
        dd = dict(row[1])

        for column, slot in list(map_column2slot.items()):
            val = dd[column]
            if ',' in val:
                val = val.split(',')
            if slot in combination:
                raise ValueError("duplicate slots: {}".format(slot))
            combination[slot] = val

        combination["name"] = "-".join([re.sub(" ", "_", dd[x]) for x in label_columns])
        combinations.append(combination)

    ret_config = {PLACEHOLDER_KEY: combinations}
    for shared_value in shared_values:
        ret_config[shared_value] = config[shared_value]

    return ret_config


def build_combinations_from_regex(config):

    slots = {}
    to_remove = set()
    for regex_key, regex_pattern in config.items():
        short_key = regex_key[:-len("_regex")]
        if regex_key.endswith("_regex") and short_key in config:
            rx = re.compile(regex_pattern)
            pairs = {}
            for filename in config[short_key]:
                try:
                    mx = rx.search(filename)
                except AttributeError as ex:
                    raise ValueError("file {} does not match '{}'".format(
                        filename, regex_pattern))
                if len(rx.groupindex) > 0:
                    key = tuple(mx.groupdict().items())
                else:
                    # naive key: position of regexs matches
                    key = (("key", "_".join(rx.search(filename).groups())), )
                pairs[key] = filename
            slots[short_key] = pairs
            to_remove.add(short_key)
            to_remove.add(regex_key)

    group_keys = collections.defaultdict(set)
    for slot_key, slot in slots.items():
        for slotitems, _ in slot.items():
            for group_key, group_value in slotitems:
                group_keys[group_key].add(group_value)

    combinations = []                
    for group_key_combination in itertools.product(*group_keys.values()):
        param_dict = dict((zip(group_keys.keys(), group_key_combination)))
        combination = {}
        name = {}
        for slot_key, slot in slots.items():
            for groups, filename in slot.items():
                for group_key, group_value in groups:
                    if group_key not in param_dict:
                        continue
                    if group_value != param_dict[group_key]:
                        break
                    name[group_key] = group_value
                else:
                    combination[slot_key] = filename

        # ignore combos where not all slots could be filled
        if len(combination) != len(slots):
            continue

        combination["name"] = "_".join(name.values())
        combinations.append(combination)

    ret_config = {PLACEHOLDER_KEY: combinations}
    for k, v in config.items():
        if k not in to_remove:
            ret_config[k] = v
    return ret_config


def build_combinations_from_config(config):

    # add multiplicity of input files
    try:
        variable = [(k, v) for k, v in list(config.items())
                    if isinstance(v, list)]
    except AttributeError:
        raise ValueError(
            "issue with configuration for option '{}'. "
            "possibly due to supplying option for tool directly and "
            "not using 'options'".format(config))

    variable = [x for x in variable if x[0] not in RESERVED_WORDS]

    combinations = []
    if variable:
        constant = [(k, v) for k, v in list(config.items())
                    if not isinstance(v, list)]
        levels = [x[0] for x in variable]
        values = [merge_shared_values(x[1]) for x in variable]
        for combination in itertools.product(*values):
            d = dict(constant + list(zip(levels, combination)))
            combinations.append(d)
    else:
        combinations.append(config)
    return combinations


def build_combinations(config):
    """build combinations of configuration parameters

    Return all possible combinations between configuration
    values. 

    There are two types of combinatorics that are applied::

       option1:
         - value1
         - value2
       option2:
         - valueA
         - valueB

    will combine into::

       - option1/value1 x option2/valueA
       - option1/value1 x option2/valueB
       ...

    Values can be grouped on lower levels for those tools expecting
    multiple input files of the same type, for a collection of samples
    to process::

       option1:
          - group1:
             - value1
             - value2
          - group2:
             - value3
             - value4
       option2:
         - valueA
         - valueB

    Will result in::
       - option1/[value1, value2] x option2/valueA
       - option1/[value1, value2] x option2/valueB
       ...

    For tools requiring multiple input files, such as a group of
    samples and a reference sequence, use the following syntax::

      test1:
         bam:
            - value1
            - value2
         reference: value10
      test2:
         bam:
            - value3
            - value4
         reference: value11

      groupby: label

    This translates into:
       - bam/[value1/values2] x reference/value10
       - bam/[value3/values4] x reference/value11

    Note the ``groupby`` variable indicationg that options
    should be grouped by the top-level (label).

    Configurations values taking multiple values are
    identified by lists, for example::

    Arg:
        config(dict) : Configuration directory

    Returns:
        list : List of dictionaries
    """

    if not config:
        return [{}]

    groupby = "option"
    if "groupby" in config:
        groupby = config["groupby"].strip()
        if groupby not in ("label", "option", "file", "regex"):
            raise ValueError(
                "unknown groupby option '{}', "
                "expected one of {}".format(
                    groupby, str(("label", "option", "file", "regex"))))
        del config["groupby"]

    if groupby == "file":
        config = build_combinations_from_design_file(config)
        groupby = "option"

    if groupby == "regex":
        config = build_combinations_from_regex(config)
        groupby = "option"

    if groupby == "option":
        combinations = build_combinations_from_config(config)
    elif groupby == "label":
        combinations = []
        for k, v in list(config.items()):
            assert isinstance(v, dict)
            combinations.append(v)

    for c in combinations:
        if PLACEHOLDER_KEY in c:
            c.update(c[PLACEHOLDER_KEY])
            del c[PLACEHOLDER_KEY]

    return sorted(combinations, key=lambda x: sorted([k, str(v)] for k, v in x.items()))
