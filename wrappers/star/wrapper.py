__author__ = "Johannes Köster"
__copyright__ = "Copyright 2016, Johannes Köster"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

assert (len(snakemake.input) == 1 or len(snakemake.input) ==2),  "input must be either 1 SE read or 2 PE reads"
fq1 = snakemake.input[0]
assert fq1 is not None, "input-> fq1 is a required input parameter"
fq1 = [snakemake.input[0]] if isinstance(snakemake.input[0], str) else snakemake.input[0]
input_str_fq1 = ",".join(fq1)
input_str = input_str_fq1
if len(snakemake.input) == 2:
    fq2 = [snakemake.input[1]] if isinstance(snakemake.input[1], str) else snakemake.input[1]
    assert len(fq1) == len(fq2), "input-> equal number of files required for fq1 and fq2"
    input_str_fq2 = ",".join(fq2) if fq2 is not None else ""
    input_str =  " ".join([input_str_fq1, input_str_fq2])


if fq1[0].endswith(".gz"):
    readcmd = "--readFilesCommand zcat"
else:
    readcmd = ""

outprefix = snakemake.output[0].split(".")[0]

# updated output to be consistent with current BBC workflow
shell("STAR "
    "{extra} "
    "--runThreadN {snakemake.threads} "
    "--genomeDir {snakemake.params.index} "
    "--readFilesIn {input_str} "
    "--twopassMode Basic "
    "{readcmd} "
    "--outSAMtype BAM SortedByCoordinate "
    "--outFileNamePrefix {outprefix}. "
    "--outStd Log "
    "{log}")
