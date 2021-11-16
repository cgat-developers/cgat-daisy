"""toolkit of a few convenience functions.

Needs to be replaced with GenomicsPLC functions.
"""

import subprocess
import contextlib
import json
import shutil
import hashlib
import signal
import collections
import logging as L
import re
import cgatcore.iotools as IOTools
import cgatcore.experiment as E
import cgatcore.pipeline as P


def run(command, return_stderr=False, **kwargs):
    """drop-in helper function to execute command line statement.
    """
    L.debug("running command: {}".format(command))

    if return_stderr:
        # expect that process fails
        p = subprocess.Popen(command, shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             **kwargs)
        stdout, stderr = p.communicate()
        return stderr
    else:
        return subprocess.check_output(command, shell=True,
                                       **kwargs)


def hash(v, encoding="utf-8"):
    """return hash value for string v"""
    return hashlib.md5(v.encode(encoding)).hexdigest()[:6]


def read_data(fn, prefix=None):
    with open(fn) as inf:
        data = json.loads(inf.read())

    if prefix is not None:
        data = dict([(prefix + key, value) for key, value in list(data.items())])

    return data


def touch(filename, times=None):
    '''update/create a sentinel file.
    Compressed files (ending in .gz) are created
    as empty 'gzip' files, i.e., with a header.
    '''
    IOTools.touch_file(filename, times=times)


def as_namedtuple(d):
    """return dictionary as a named tuple.
    """
    return collections.namedtuple('GenericDict', list(d.keys()))(**d)


def update_namedtuple(t, **kwargs):
    """add additional fields to a named tuple

    returns a new tuple.
    """
    a = t._asdict()
    a.update(kwargs)
    return collections.namedtuple('GenericDict', list(a.keys()))(**a)


def parse_region_string(s):
    """parse a genomic region string.

    Returns tuple of contig, start, end. Missing values are None.
    """
    if not s:
        return None, None, None
    if ":" in s:
        contig, coords = s.split(":")
        if "-" in coords:
            start, end = list(map(int, coords.split("-")))
        else:
            start = int(start)
            end = None
        return contig, start, end
    else:
        return s, None, None


def sort_by_chromosome(strings):

    prog = re.compile(".*(chr[\dXYMT]+).*")

    def extract_chrom_string(string):
        match = prog.match(string)
        if match:
            return match.group(1)
        else:
            return ''

    max_chrom_number = 50
    map_dict2order = dict(
        ["chr{}".format(y), x] for x, y in enumerate(range(max_chrom_number))
    )

    map_dict2order["chrX"] = max_chrom_number + 1
    map_dict2order["chrY"] = max_chrom_number + 2
    map_dict2order["chrMT"] = max_chrom_number + 2

    return sorted(strings, key=lambda x: map_dict2order.get(
        extract_chrom_string(x), 1000))


def redirect2mounts(config,
                    mountpoint=None,
                    debug=None,
                    mount_write=False,
                    substitute_only=False,
                    always_mount=False):
    """redirect filenames in dictionary config to a mount-point.

    Mount points in the config are indicated by the `arv=` prefix. If
    no option in config requires mounting, no mounting will be done and
    the method returns None.

    :param config: dictionary with config values. Will be modified in-place.
    :param mountpoint: if given, paths will be substituted by mountpoint. If None,
        a new mountpoint will be created.
    :param debug: if given, mount in debug mode and save log to filename.
    :param mount_write: if True, mount in --read-write mode.
    :param substitute_only: if True, only perform substitution, do not mount anything
        even if mountpoint is None.
    :param always_mount: if True, always mount, no matter if arv= prefix is present.

    :return: the mountpoint

    """
    arvados_options = ["--disable-event-listening"]
    if debug:
        arvados_options.append(" --debug --logfile={}".format(debug))

    if mount_write:
        arvados_options.append("--read-write")
        arvados_options = " ".join(arvados_options)
        if not mountpoint:
            mountpoint = P.get_temp_dir() + "/"
            E.info("redirect2mounts: mounting arvados at {} with --read-write".format(
                mountpoint))
            E.run("arv-mount {} {}".format(
                arvados_options, mountpoint))
            E.info("redirect2mounts: arvados mounted at {} with --read-write".format(
                mountpoint))
    else:
        arvados_options.append("--read-only")
        if always_mount:
            mountpoint = P.get_temp_dir() + "/"
            do_mount = True
        else:
            do_mount = False

        for d, key, value in IOTools.nested_iter(config):
            if isinstance(value, str):
                if "arv=" in value:
                    if substitute_only and mountpoint is None:
                        continue
                    if not mountpoint:
                        mountpoint = P.get_temp_dir() + "/"
                        do_mount = True
                    d[key] = re.sub("arv=", mountpoint, value)

        if do_mount:
            raise NotImplementedError("arvados support disabled")
            # if not arvados.have_arvados():
            #     raise ValueError(
            #         "config file requires arvados access, but arvados not available")
            arvados_options = " ".join(arvados_options)
            E.debug("redirect2mounts: mounting arvados at {} with options {}".format(
                mountpoint, arvados_options))
            E.run("arv-mount {} {}".format(
                arvados_options, mountpoint))
            E.debug("redirect2mounts: arvados mounted at {}".format(
                mountpoint))

    return mountpoint


@contextlib.contextmanager
def arvados_enabled(*args, **kwargs):
    """accepts the same options as redirect2mounts"""

    def arvados_cleanup(signal=None, frame=None):

        mountpoint = P.PARAMS.get("mount_point", None)
        if mountpoint:
            E.debug("unmounting arvados at {}".format(mountpoint))
            E.run("fusermount -u {}".format(mountpoint))
            E.debug("unmounted arvados at {}".format(mountpoint))
            try:
                shutil.rmtree(mountpoint)
            except OSError:
                # ignore errors - arvados issue (read-only file system)?
                E.warn("failure while removing mountpoint {} - ignored".format(mountpoint))
            P.PARAMS["mount_point"] = None

    signal.signal(signal.SIGINT, arvados_cleanup)
    args = list(args)
    if args:
        params = args.pop()
    else:
        params = kwargs.get("config", P.PARAMS)

    mountpoint = redirect2mounts(params, *args, **kwargs)

    if mountpoint is not None:
        params["mount_point"] = mountpoint
        P.PARAMS["mount_point"] = mountpoint

    try:
        yield params
    finally:
        E.debug("cleaning up arvados")
        arvados_cleanup()
