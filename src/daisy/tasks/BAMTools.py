import glob
import os
import pipes
import math
import re
import shutil
import tempfile
from .ToolRunner import ToolRunner
from .CollateRunner import CollateRunner
from .Runner import resolve_argument
from .BAMMetrics import BAMPreprocessor
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


def sra_peek(sra, outdir=None):
    """return the full file names for all files which will be extracted

    Parameters
    ----------

    outdir : path
        perform extraction in outdir. If outdir is None, the extraction
        will take place in a temporary directory, which will be deleted
        afterwards.

    Returns
    -------
    files : list
        A list of fastq formatted files that are contained in the archive.
    format : string
        The quality score format in the :term:`fastq` formatted files.

    """

    if outdir is None:
        workdir = tempfile.mkdtemp()
    else:
        workdir = outdir

    # --split-files creates files called prefix_#.fastq.gz,
    # where # is the read number.
    # If file cotains paired end data:
    # output = prefix_1.fastq.gz, prefix_2.fastq.gz
    #    *special case: unpaired reads in a paired end --> prefix.fastq.gz
    #    *special case: if paired reads are stored in a single read,
    #                   fastq-dump will split. There might be a joining
    #                   sequence. The output would thus be:
    #                   prefix_1.fastq.gz, prefix_2.fastq.gz, prefix_3.fastq.gz
    #                   You want files 1 and 3.

    E.run("""fastq-dump --split-files --gzip -X 1000
                 --outdir %(workdir)s %(sra)s""" % locals())
    f = sorted(glob.glob(os.path.join(workdir, "*.fastq.gz")))
    ff = [os.path.basename(x) for x in f]

    if len(f) == 1:
        # sra file contains one read: output = prefix.fastq.gz
        pass

    elif len(f) == 2:
        # sra file contains read pairs:
        # output = prefix_1.fastq.gz, prefix_2.fastq.gz
        assert ff[0].endswith(
            "_1.fastq.gz") and ff[1].endswith("_2.fastq.gz")

    elif len(f) == 3:
        if ff[2].endswith("_3.fastq.gz"):
            f = glob.glob(os.path.join(workdir, "*_[13].fastq.gz"))
        else:
            f = glob.glob(os.path.join(workdir, "*_[13].fastq.gz"))

    if outdir is None:
        shutil.rmtree(workdir)

    return [os.path.basename(x) for x in f]


def build_readgroup_string(outfile, params):

    if params.readgroup_id_regex is None:
        readgroup_id = IOTools.snip(os.path.basename(outfile), ".bam")
    else:
        try:
            readgroup_id = "-".join(re.search(
                params.readgroup_id_regex,
                outfile).groups())
        except AttributeError as ex:
            raise AttributeError("regular expression {} does not match {}".format(
                params.readgroup_id_regex, outfile))

    if params.readgroup_sample_regex is None:
        readgroup_sample = readgroup_id
    else:
        try:
            readgroup_sample = "-".join(re.search(
                params.readgroup_sample_regex,
                outfile).groups())
        except AttributeError as ex:
            raise AttributeError("regular expression {} does not match {}".format(
                params.readgroup_sample_regex, outfile))

    readgroup_string = "@RG\tID:{}\tSM:{}".format(
        readgroup_id, readgroup_sample)

    if params.readgroup_header:
        readgroup_string += "\t{}".format(params.readgroup_header)

    return readgroup_string, readgroup_id, readgroup_sample


class run_collate_samtools_merge(CollateRunner):
    path = "samtools"
    name = "samtools_merge"

    set_readgroup = False
    readgroup_id_regex = None
    readgroup_sample_regex = None
    readgroup_header = ""

    # options for samtools view
    view_options = ""

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, infiles, outfile, params):

        files = " ".join(infiles)
        job_threads = 1
        if "--threads" in params.options:
            job_threads = int(re.search("--threads[= ]\s*(\d+)",
                                        params.options).groups()[0])

        if "--threads" in params.view_options:
            job_threads += int(re.search("--threads[= ]\s*(\d+)",
                                         params.view_options).groups()[0])

        if params.set_readgroup or params.readgroup_id_regex is not None:
            readgroup_string, readgroup_id, readgroup_sample = build_readgroup_string(
                outfile, params)

            with open(outfile + ".header.sam", "w") as outf:
                outf.write(readgroup_string + "\n")

            retval = P.run(
                "{params.path} merge "
                "{params.options} "
                "-f "
                "-h {outfile}.header.sam "
                "-r "
                "- "
                "{files} "
                "2> {outfile}.log "
                "| samtools view -h - "
                "| perl -p -e 's/^.*\\n// if (/^\@RG/ && !/{readgroup_id}/); "
                "   s/RG:Z:\S+/RG:Z:{readgroup_id}/' "
                "| samtools view -bS "
                "{params.view_options} "
                "- "
                "> {outfile}; "
                "samtools index {outfile} 2> {outfile}.index.log"
                .format(**locals()), job_threads=job_threads)

        else:
            retval = P.run(
                "{params.path} merge "
                "-f "
                "{params.options} "
                "{outfile} "
                "{files} "
                "2> {outfile}.log; "
                "samtools index {outfile} 2> {outfile}.index.log"
                .format(**locals()), job_threads=job_threads)
        return retval


class BWARunner(ToolRunner):
    output = "result.bam"
    path = "bwa"

    def get_version(self):
        help_text = E.run("{self.path}".format(**locals()),
                          return_stderr=True).strip()
        return re.search("Version: (\S+)", help_text).groups()[0]


class run_tool_bwa_mem(BWARunner):
    name = "bwa_mem"
    expected = ["reference_fasta", "fastq"]

    set_readgroup = False
    readgroup_id_regex = None
    readgroup_sample_regex = None
    readgroup_header = ""

    def run(self, outfile, params):

        if "-t" in params.options:
            job_threads = int(re.search("-t\s*(\d+)",
                                        params.options).groups()[0])
        else:
            job_threads = 1

        # BWA requires at least 6Gb of memory, but is also correlated
        # with the number of threads, so use 5Gb + 1Gb per thread
        job_memory = "{}G".format(5.0 + 1.0 * job_threads)

        fastq = resolve_argument(params.fastq, ",")
        fastq = '"{}"'.format('" "'.join(fastq.split(",")))

        tmpdir = P.get_temp_filename(clear=True)

        if params.set_readgroup or params.readgroup_id_regex is not None:
            readgroup_string, readgroup_id, readgroup_sample = build_readgroup_string(
                outfile, params)

            # pipes.quote needs to shlex.quote in py3
            readgroup_option = "-R {}".format(pipes.quote(readgroup_string))
            # add additional level of quoting:
            readgroup_option = re.sub("\\t", "\\\\t", readgroup_option)
        else:
            readgroup_option = ""

        return P.run(
            "mkdir {tmpdir}; "
            "{self.path} mem "
            "{readgroup_option} "
            "{params.options} "
            "{params.reference_fasta} "
            "{fastq} "
            "2> {outfile}.log "
            "| samtools view -bu /dev/stdin "
            "2> {outfile}.view.log "
            "| samtools sort --threads {job_threads} -T {tmpdir} -O bam /dev/stdin "
            "2> {outfile}.sort.log "
            "> {outfile}; "
            "samtools index {outfile} >& {outfile}.index.log; "
            "rm -rf {tmpdir}".format(**locals()),
            **params._asdict())


class run_tool_bwa_mem_sra(BWARunner):
    """extract fastq sequences from SRA files and map with bwa.

    extract_to_temp: bool
       if set to False, the processing will happen inside the
       current working directory. The default is to use a temporary
       directory. For protected data, use extract_to_temp = False.

    """

    name = "bwa_mem_sra"
    output = "result.cram"

    path_fastqdump = "fastq-dump"
    expected = ["reference_fasta", "sra"]

    set_readgroup = False
    readgroup_id_regex = None
    readgroup_sample_regex = None
    readgroup_header = ""

    cram_fasta = None

    extract_to_temp = True

    def run(self, outfile, params):

        min_job_memory = 3
        if "-t" in params.options:
            job_threads = int(re.search("-t\s*(\d+)",
                                        params.options).groups()[0])
        else:
            job_threads = 1

        job_memory = "{}G".format(
            float(min_job_memory + 1.0 * job_threads) / job_threads)

        cram_fasta = params.cram_fasta
        if params.cram_fasta is None:
            cram_fasta = params.reference_fasta

        if params.set_readgroup or params.readgroup_id_regex is not None:
            readgroup_string, readgroup_id, readgroup_sample = build_readgroup_string(
                outfile, params)

            # pipes.quote needs to shlex.quote in py3
            readgroup_option = "-R {}".format(pipes.quote(readgroup_string))
            # add additional level of quoting:
            readgroup_option = re.sub("\\t", "\\\\t", readgroup_option)
        else:
            readgroup_option = ""

        fastq = " ".join(sra_peek(params.sra))
        outfile = os.path.abspath(outfile)

        if params.extract_to_temp:
            tmpdir = P.get_temp_filename(clear=True)
            tmpdir_pre = "mkdir {};".format(tmpdir)
            tmpdir_post = "rm -rf {}".format(tmpdir)
        else:
            tmpdir = os.path.dirname(outfile)
            tmpdir_pre = ""
            tmpdir_post = ""

        # AH: fastq-dump hangs with arv mounts, thus try copying first
        if not IOTools.is_local(params.sra):
            E.warn("copying file {} to temporary directory".format(params.sra))
            temp_sra = os.path.join(
                tmpdir, os.path.basename(params.sra))
            fastq_dump = (
                "cp {params.sra}* {tmpdir}; "
                "fastq-dump --split-files --gzip {temp_sra} >& {outfile}.dump.log ".format(
                    **locals()))
            tmpdir_post = "rm -f {}*; {}".format(
                temp_sra, tmpdir_post)
        else:
            fastq_dump = (
                "fastq-dump --split-files --gzip {params.sra} >& {outfile}.dump.log "
            )

        return P.run(
            "{tmpdir_pre} "
            "cd {tmpdir}; "
            "{fastq_dump}; "
            "{self.path} mem -v 3 "
            "{readgroup_option} "
            "{params.options} "
            "{params.reference_fasta} "
            "{fastq} "
            "2> {outfile}.map.log "
            "| samtools view -O cram --reference {params.cram_fasta} /dev/stdin "
            "2> {outfile}.view.log "
            "| samtools sort -T {tmpdir} -O cram /dev/stdin "
            "2> {outfile}.sort.log "
            "> {outfile}; "
            "samtools index {outfile} >& {outfile}.index.log; "
            "{tmpdir_post}".format(**locals()))


class PicardAddReadGroup():
    path_picard = "/data/install/Free/picard-tools-1.140/picard.jar"
    sample = "unknown"
    platform = "illumina"
    library = "unknown"

    def build_picard_statement(self, infile, outfile, params):

        if params.set_readgroup or params.readgroup_id_regex is not None:
            readgroup_string, readgroup_id, readgroup_sample = build_readgroup_string(
                outfile, params)

            readgroup_options = (
                "RGID=1 "
                "RGLB={params.library} "
                "RGPL={params.platform} "
                "RGPU=unknown "
                "RGSM={readgroup_sample} ".format(**locals()))
        else:
            readgroup_options = (
                "RGID=1 "
                "RGLB={params.library} "
                "RGPL={params.platform} "
                "RGPU=unknown "
                "RGSM={params.sample} ".format(**locals()))

        statement = (
            "java -Xmx8000m -jar {params.path_picard} "
            "AddOrReplaceReadGroups "
            "INPUT={infile} "
            "OUTPUT={outfile} "
            "VALIDATION_STRINGENCY=LENIENT "
            "{readgroup_options} "
            ">& {outfile}.picard.log; "
            "samtools index {outfile}".format(**locals()))
        return statement


class run_tool_isaac_align(BWARunner, PicardAddReadGroup):
    name = "isaac_align"
    expected = ["reference_fasta", "fastq"]

    path = "isaac-align"

    bam_gzip_level = 6

    set_readgroup = False
    readgroup_id_regex = None
    readgroup_sample_regex = None
    readgroup_header = ""

    def get_version(self):
        help_text = E.run("{self.path} --version".format(**locals()),
                          return_stdout=True).strip()
        return help_text

    def run(self, outfile, params):

        local_options = []
        outfile = os.path.abspath(outfile)
        outdir = os.path.dirname(outfile)

        # assumption is that index is called xyz.fa without the .fa.
        reference_fasta = IOTools.snip(params.reference_fasta, ".fa", ".fasta")
        if not os.path.exists(reference_fasta):
            raise ValueError("input reference {} does not exist".format(reference_fasta))

        if "--jobs" in params.options or "-j" in params.options:
            job_threads = int(re.search("(--jobs|-j)\s*(\d+)",
                                        params.options).groups()[1])
        else:
            job_threads = 8

        if "--memory-limit" in params.options or "-m" in params.options:
            job_memory_gb = int(re.search("(--memory-limit|-m)\s*(\d+)",
                                          params.options).groups()[1])
        else:
            job_memory_gb = 60
            local_options.append("--memory-limit {}".format(job_memory_gb))

        if job_memory_gb < 60:
            E.warn("isaac-align likely to require at least 60Gb of memory, {}G requested".format(
                job_memory_gb))

        job_memory = "{}G".format(float(job_memory_gb) / job_threads)

        fastq_dir = os.path.join(outdir, "input_fastq")
        if not os.path.exists(fastq_dir):
            os.makedirs(fastq_dir)

        if len(params.fastq) == 2:
            if not os.path.exists(os.path.join(fastq_dir, "lane1_read1.fastq.gz")):
                os.symlink(os.path.abspath(params.fastq[0]), os.path.join(fastq_dir, "lane1_read1.fastq.gz"))
            if not os.path.exists(os.path.join(fastq_dir, "lane1_read2.fastq.gz")):
                os.symlink(os.path.abspath(params.fastq[1]), os.path.join(fastq_dir, "lane1_read2.fastq.gz"))
        else:
            raise NotImplementedError("expected 2 fastq files, received only {}".format(len(params.fastq)))

        intermediate_bam = os.path.join(outdir,
                                        "Aligned",
                                        "Projects",
                                        "default",
                                        "default",
                                        "sorted.bam")

        # picard statement to set readgroup
        picard_statement = self.build_picard_statement(
            intermediate_bam,
            outfile,
            params)

        tmpdir = os.path.join(outdir, "TEMP")

        local_options = " ".join(local_options)
        # isaac generates output files in working directory, so do a cd and make
        # sure that absolute path names are used elsewhere.
        statement = (
            "cd {outdir}; "
            "{self.path} "
            "--reference-genome {reference_fasta}/sorted-reference.xml "
            "--base-calls {fastq_dir} "
            "--base-calls-format fastq-gz "
            "--temp-directory {tmpdir} "
            "--cleanup-intermediary 1 "
            "--bam-gzip-level {params.bam_gzip_level} "
            "{params.options} "
            "{local_options} "
            ">& {outfile}.isaac.log; "
            "{picard_statement}; "
            "rm -rf {tmpdir} "
            .format(**locals()))

        return P.run(statement)


class run_tool_graphmap(BWARunner, PicardAddReadGroup):
    name = "graphmap"
    expected = ["reference_fasta", "fastq"]
    path = "graphmap"

    def get_version(self):
        help_text = E.run("{self.path} -h".format(**locals()),
                          return_stdout=True).strip()
        return re.search("Version: (\S+)", help_text).groups()[0]

    def run(self, outfile, params):

        if "-t" in params.options:
            job_threads = int(re.search("-t\s*(\d+)",
                                        params.options).groups()[0])

        job_memory = "32G"

        fastq = resolve_argument(params.fastq, " ")

        tmpdir = P.get_temp_filename(clear=True)

        return P.run(
            "mkdir {tmpdir}; "
            "{params.path} "
            "{params.options} "
            "-r {params.reference_fasta} "
            "-d {fastq} "
            "-o {tmpdir}/result.sam "
            ">& {outfile}.log; "
            "samtools view -bS {tmpdir}/result.sam "
            "| samtools sort -o {tmpdir}/sorted.bam -; "
            "java -Xmx8000m -jar {params.path_picard} "
            "AddOrReplaceReadGroups "
            "INPUT={tmpdir}/sorted.bam "
            "OUTPUT={outfile} "
            "VALIDATION_STRINGENCY=LENIENT "
            "RGID=1 "
            "RGLB={params.library} "
            "RGPL={params.platform} "
            "RGPU=unknown "
            "RGSM={params.sample} "
            ">& {outfile}.picard.log; "
            "samtools index {outfile} "
            ">& {outfile}.index.log; "
            "rm -rf {tmpdir}".format(**locals()))


class run_tool_bbmap(BWARunner, PicardAddReadGroup):
    name = "bbmap"
    expected = ["reference_fasta", "fastq"]
    path = "bbmap"

    def get_version(self):
        help_text = E.run("{self.path} -version".format(**locals()),
                          return_stderr=True).strip()
        if help_text and "not found" not in help_text:
            return re.search(r"BBMap version (\S+)", help_text).groups()[0]
        else:
            raise ValueError("bbmap not found at/as {}: {}".format(
                self.path, help_text))

    def run(self, outfile, params):

        # the default is auto so use ten threads.
        if "threads" in params.options:
            if "job_threads=auto" in params.options:
                raise ValueError(
                    "please specify the number of threads "
                    "to use explicitely")
            else:
                job_threads = int(re.search("threads=(\d+)",
                                            params.options).groups()[0])
        else:
            raise ValueError("please specify the number of threads to use")

        job_memory = "32G"

        fastq = resolve_argument(params.fastq, " ")

        tmpdir = P.get_temp_filename(clear=True)

        return P.run(
            "mkdir {tmpdir}; "
            "zcat {fastq} "
            "| cut -c -5999 "
            "| gzip > {tmpdir}/in.fastq.gz; "
            "{params.path} "
            "{params.options} "
            "in={tmpdir}/in.fastq.gz "
            "ref={params.reference_fasta} "
            "out={tmpdir}/result.bam "
            ">& {outfile}.log; "
            "samtools sort -o {tmpdir}/sorted.bam {tmpdir}/result.bam; "
            "java -Xmx8000m -jar {params.path_picard} "
            "AddOrReplaceReadGroups "
            "INPUT={tmpdir}/sorted.bam "
            "OUTPUT={outfile} "
            "VALIDATION_STRINGENCY=LENIENT "
            "RGID=1 "
            "RGLB={params.library} "
            "RGPL={params.platform} "
            "RGPU=unknown "
            "RGSM={params.sample} "
            ">& {outfile}.picard.log; "
            "samtools index {outfile} "
            ">& {outfile}.index.log; "
            "rm -rf {tmpdir}".format(**locals()))


class BowtieRunner(ToolRunner):
    output = "result.bam"
    path = "bowtie"

    def get_version(self):
        help_text = E.run("{self.path} --version 2> /dev/null".format(**locals()),
                          return_stdout=True).strip()
        return re.search("version (\S+)", help_text).groups()[0]


class run_tool_bowtie2(BowtieRunner):
    name = "bowtie2"
    path = "bowtie2"

    expected = ["reference_fasta", "fastq"]

    set_readgroup = False
    readgroup_id_regex = None
    readgroup_sample_regex = None
    readgroup_header = ""

    def run(self, outfile, params):

        if "--threads" in params.options or "-t " in params.options:
            job_threads = int(re.search("(-t|--threads)\s*(\d+)",
                                        params.options).groups()[1])

        fastq = resolve_argument(params.fastq, ",").split(",")
        if len(fastq) == 1:
            fastq = '-U "{}"'.format(fastq)
        else:
            fastq = '-1 "{}" -2 "{}"'.format(*fastq)

        tmpdir = P.get_temp_filename(clear=True)

        if "index" in params._fields:
            index = params.index
        else:
            index = params.reference_fasta

        if params.set_readgroup or params.readgroup_id_regex is not None:
            readgroup_string, readgroup_id, readgroup_sample = build_readgroup_string(
                outfile, params)

            # pipes.quote needs to shlex.quote in py3
            readgroup_option = "--rg-id {}".format(readgroup_id)

            # add additional level of quoting and remove "ID:{}"
            readgroup_string = re.sub("@RG\tID:\S+\t", "", readgroup_string)
            readgroup_string = " ".join(["--rg {}".format(x)
                                         for x in readgroup_string.split("\t")])
        else:
            readgroup_option = ""
            readgroup_string = ""

        return P.run(
            "mkdir {tmpdir}; "
            "{self.path} "
            "{readgroup_option} "
            "{readgroup_string} "
            "{params.options} "
            "-x {index} "
            "{fastq} "
            "2> {outfile}.log "
            "| samtools view -b /dev/stdin "
            "2> {outfile}.view.log "
            "| samtools sort -T {tmpdir} -O bam /dev/stdin "
            "2> {outfile}.sort.log "
            "> {outfile}; "
            "samtools index {outfile}; "
            "rm -rf {tmpdir}".format(**locals()),
            **params._asdict())


class run_tool_import_bam(ToolRunner, BAMPreprocessor):

    name = "import_bam"

    expected = ["bam"]
    output = "result.bam"

    def get_version(self):
        return "builtin"

    def run(self, outfile, params):

        return self.run_with_preprocessing(
            params.bam,
            outfile,
            params,
            "cp {params.bam} {outfile}; "
            "cp {params.bam}.bai {outfile}.bai".format(**locals()),
            **params._asdict())


class run_tool_picard_markduplicates(ToolRunner):
    expected = ["reference_fasta", "bam"]
    output = "result.bam"
    path = "/data/install/Free/picard-tools-1.140/picard.jar"
    name = "picard_markduplicates"

    def get_version(self):
        return E.run(
            "java -jar {self.path} SortVcf --version".format(**locals()),
            return_stderr=True).strip()

    def run(self, outfile, params):

        bam = resolve_argument(params.bam)

        # rename index from x.bai to x.bam.bai
        outprefix = IOTools.snip(outfile, ".bam", ".cram")

        statement = ("java -Xmx8000m -jar {params.path} "
                     "MarkDuplicates "
                     "INPUT={bam} "
                     "TMP_DIR=%(tmpdir)s "
                     "CREATE_INDEX=TRUE "
                     "REFERENCE_SEQUENCE={params.reference_fasta} "
                     "METRICS_FILE={outfile}.metrics "
                     "{params.options} "
                     "OUTPUT={outfile} "
                     ">& {outfile}.log; "
                     "mv {outprefix}.bai {outfile}.bai".format(**locals()))

        # 12G is required for java overhead
        return P.run(statement, job_memory="12G")
