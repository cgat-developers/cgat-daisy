"""Tool and metric library
"""
import inspect
import os
import glob
import importlib
import importlib.util
# bioconda: ruamel.yml
# defaults: ruamel_yaml
try:
    import ruamel_yaml as yaml
    from ruamel_yaml import RoundTripLoader
except ImportError:
    from ruamel import yaml
    from ruamel.yaml import RoundTripLoader

import cgatcore.iotools as IOTools
import daisy.toolkit
import cgatcore.experiment as E


@E.cached_function
def get_default_params():
    """return default parameters for tools/metrics.

    Could be refactored to read defaults from a user specified file.
    The current implementation takes the one located within the
    repository.
    """

    defaults = os.path.join(os.path.dirname(__file__), "defaults.yml")
    if os.path.exists(defaults):
        with IOTools.open_file(defaults) as inf:
            result = yaml.load(inf, Loader=RoundTripLoader)
    else:
        result = {}
    return result


def redirect_defaults2mountpoint(mountpoint):
    """redirect paths in the defaults to mountpoint.

    Redirection only happens if the config session indicates that
    arvados is present.
    """
    params = get_default_params()
    mountpoint = daisy.toolkit.redirect2mounts(params,
                                               mountpoint,
                                               substitute_only=True)
    return mountpoint


def get_task_version(taskf):
    """return a version for a taskf.

    This is wrapper around taskf.get_version() that will not raise
    an exception if a tool is not installed but will return None
    instead.
    """
    try:
        version = taskf().get_version()
    except Exception as ex:
        msg = str(ex)
        if msg.startswith("no version defined for"):
            raise ValueError(msg)
        elif "non-zero exit status" in msg:
            # if tool can not be found, that is ok
            return None
        elif "not found" in msg:
            return None
        else:
            raise ValueError("version check failed with: {}".format(msg))

    return version

module_dirs = [os.path.join(os.path.dirname(__file__))]
module_dirs.extend([x.strip() for x in
                    os.environ.get("DAISY_TASKLIBRARY", "").split(",")
                    if x.strip()])

modules = []
for idx, root in enumerate(module_dirs):
    for module in glob.glob(os.path.join(root, "*.py")):
        if "flycheck" in module:
            continue
        if module.endswith("__init__.py"):
            continue
        module_name = IOTools.snip(os.path.basename(module))
        if idx == 0:
            modules.append(importlib.import_module(
                "daisy.tasks.{}".format(module_name)))
        else:
            spec = importlib.util.spec_from_file_location(
                "daisy.UserLibrary.{}".format(module_name),
                module)
            foo = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(foo)
            modules.append(foo)

# TODO: use derivation instead of name prefix
map_tool_to_runner = dict()
map_metric_to_runner = dict()
map_collate_to_runner = dict()
map_split_to_runner = dict()

for module in modules:
    map_tool_to_runner.update(dict([
        (y.name, y) for x, y in inspect.getmembers(
            module,
            lambda a: inspect.isclass(a)) if x.startswith("run_tool")]))

    map_metric_to_runner.update(dict([
        (y.name, y) for x, y in inspect.getmembers(
            module,
            lambda a: inspect.isclass(a)) if x.startswith("run_metric")]))

    map_collate_to_runner.update(dict([
        (y.name, y) for x, y in inspect.getmembers(
            module,
            lambda a: inspect.isclass(a)) if x.startswith("run_collate")]))

    map_split_to_runner.update(dict([
        (y.name, y) for x, y in inspect.getmembers(
            module,
            lambda a: inspect.isclass(a)) if x.startswith("run_split")]))


__all__ = ["map_tool_to_runner",
           "map_metric_to_runner",
           "map_collate_to_runner",
           "map_split_to_runner"]
