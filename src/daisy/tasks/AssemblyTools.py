import os
import re

from .ToolRunner import ToolRunner
import cgatcore.pipeline as P
import cgatcore.experiment as E


class run_tool_canu(ToolRunner):
    name = "canu"
    path = "canu"
    path_java = "java"

    expected = ["fasta"]
    output = "result.fasta"

    genome_size = "3g"
    assembly_mode = "-nanopore-raw"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stderr=True).strip()
        if help_string and "not found" not in help_string:
            return re.search("Canu (.+)", help_string).groups()[0]
        else:
            raise ValueError("canu not found at/as {}: {}".format(
                self.path, help_string))

    def run(self, outfile, params):

        path = os.environ["PATH"]
        gp = P.get_parameters_as_namedtuple()
        cluster_queue = gp.cluster["queue"]
        cluster_memory_resource = gp.cluster["memory_resource"]
        cluster_parallel_environment = gp.cluster["parallel_environment"]
        outdir = os.path.dirname(outfile)
        outname = os.path.basename(outdir)
        # -sync y forces qsub to wait until job completes before
        # continuing.

        statement = (
            "{self.path} "
            "-p canu "
            "-d {outdir} "
            "-genomeSize={params.genome_size} "
            "gridOptionsJobName={outname} "
            "java={params.path_java} "
            "gridOptions=\"-q {cluster_queue} -v PATH={path} -sync y \" "
            "gridEngineMemoryOption=\"-l {cluster_memory_resource}=MEMORY\" "
            "gridEngineThreadsOption=\"-pe {cluster_parallel_environment} THREADS\" "
            "{params.options} "
            "{params.assembly_mode} "
            "{params.fasta} "
            ">& {outfile}.log; "
            "mv {outdir}/canu.contigs.fasta {outfile}".format(**locals()))

        return P.run(statement, without_cluster=True)
