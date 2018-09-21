
import os
import sys
import glob
import inspect
import re
import importlib

from cgatcore.iotools import snip

IGNORE = ["__init__.py", "cli.py"]
SYNONYMS = {"benchmark-simple": ["run"]}


def main(argv=None):

    if argv is None:
        argv = sys.argv

    modules = []
    for module in glob.glob(os.path.join(os.path.dirname(__file__), "*.py")):
        if os.path.basename(module) in IGNORE:
            continue
        if "flycheck" in module:
            continue
        mod = "daisy.tools.{}".format(snip(os.path.basename(module)))
        modules.append(importlib.import_module(mod))

    script_dict = {}
    for module in modules:
        try:
            f = [y for (x, y) in inspect.getmembers(module) if x == "main"][0]
        except IndexError:
            continue
        name = re.sub("_", "-", module.__name__.split(".")[-1])
        script_dict[name] = f
        for synonym in SYNONYMS.get(name, []):
            script_dict[synonym] = f

    if len(argv) == 1:
        print('\n'.join(sorted(script_dict.keys())))
    else:
        command_key = argv[1]
        command_args = argv[1:]

        command = script_dict[command_key]
        try:
            return command(command_args)
        except:
            print('When running {!r}'.format(command_key))
            raise


if __name__ == '__main__':
    sys.exit(main(sys.argv[:]))
