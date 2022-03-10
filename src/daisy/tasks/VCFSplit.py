import re
import os
import pysam
import shutil
import pandas

from .SplitRunner import SplitRunner
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


class run_split_vcf_by_chromosome(SplitRunner):
    """split a VCF by chromosome using tabix.

    The file is assumed to be sorted and indexed.
    """
    name = "vcf_by_chromosome"

    output = "chromosome_{}.dir/result.vcf.gz"

    def get_version(self):
        return "builtin"

    def run(self, infile, outfiles, params):

        tbxfile = pysam.VariantFile(infile)
        statements = []

        for chrom in list(tbxfile.header.contigs):
            output_file = outfiles.format(chrom)
            output_dir = os.path.dirname(output_file)
            statements.append(
                "mkdir {output_dir}; "
                "tabix -h {infile} {chrom} | bgzip > {output_file}; "
                "tabix -p vcf {output_file} "
                .format(**locals()))

        retvals = P.run(statements)

        # clean up empty vcfs, opening empty VCF in pysam throws
        # ValueError
        for chrom in list(tbxfile.header.contigs):
            output_file = outfiles.format(chrom)
            output_dir = os.path.dirname(output_file)
            try:
                f = pysam.VariantFile(output_file)
                f.close()
            except ValueError:
                E.warn("removing empty VCF {}".format(output_file))
                shutil.rmtree(output_dir)

        tbxfile.close()
