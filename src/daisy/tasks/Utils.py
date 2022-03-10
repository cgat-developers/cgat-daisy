"""Utils

Utilities for tasks in the tasks.

"""

import os
import pysam
import cgatcore.experiment as E


@E.cached_function
def get_sequence_length_dict(fastafn):
    """return sequence/length dictionary from a
    fasta file.

    The fasta file needs to be indexed with samtools faidx.
    """
    # Temporary fix: see issue SYS-517
    if not os.path.exists(fastafn):
        E.warn("could not find file {}".format(fastafn))
        return None

    try:
        with pysam.FastaFile(fastafn) as inf:
            fastadict = dict(list(zip(
                inf.references,
                inf.lengths)))
    except IOError as ex:
        E.warn("file {} could not be opened".format(ex))
        fastadict = None

    return fastadict


def sequence_length_dicts_iterator(fastafiles):
    """iterate over sequence/length dictionary from a list of fasta files.
    """

    for fastafn in fastafiles:
        fastadict = get_sequence_length_dict(fastafn)
        if fastadict:
            yield fastafn, fastadict


def match_sequence_dictionaries(sequence_dict, fastafiles):
    """match a sequence dictionary (contig, length) against
    a collection of fasta files.

    :param sequence_dict: dictionary of contig/length pairs.
    :param fastafiles: list of :term:`fasta` formatted files. The
       fasta files need to indexed with samtools faidx.
    :return: a tuple (fastafn, diffs). Fastafn is the filename of
       the fasta file that has been matched, None if no match
       has been found. Diffs contains a list of
       discrepancies between sequence_dict and files in
       fastafiles that have been examined. If sequence_dict is empty,
       fastafn will be None and the list of diffs empty.
    """
    if not sequence_dict:
        return None, []

    fastafn = None
    diffs = []
    # match by sequence dictionary with optional length
    for fastafn, fastadict in sequence_length_dicts_iterator(fastafiles):
        E.debug("inspected {}".format(fastafn))
        contig_missing = None
        length_mismatch = None
        for reference, length in sequence_dict.items():
            if reference not in fastadict:
                contig_missing = reference
                break
            if length > 0 and length != fastadict[reference]:
                length_mismatch = reference
                break

        if not (length_mismatch or contig_missing):
            break
        else:
            diffs.append((fastafn, contig_missing, length_mismatch))

    return fastafn, diffs
