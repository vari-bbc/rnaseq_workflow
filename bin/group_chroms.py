#!/usr/bin/env python3

import os
import pandas as pd
from snakemake.io import load_configfile

config = load_configfile("config.yaml")

# make a file with col1 as contig group name and col2 as a commma separated list of contigs. This was written for use with GATK for splitting the variant calling steps by chromosome/contig. We group the unplaced contigs together since those tend to be small.
fai_file = config["ref"]["fai"]
contigs_file = "grouped_contigs.tsv"

os.system("perl -e 'print qq/name\tcontigs\n/' >" + contigs_file)
os.system("grep -Po '(^chr\S+|^\d+)' " + fai_file + " | perl -lne 'print qq/$_\t$_/' >>" + contigs_file)
os.system("grep -Pv '(^chr\S+|^\d+)' " + fai_file + " | cut -f1 | perl -npe 's/\n/,/g' | perl -lne 's/,$//; print qq/unplaced_contigs\t$_/'  >>" + contigs_file)
contig_groups = pd.read_table(contigs_file)

# check chromosomes/contigs parsed correctly.
num_contigs_fai = sum(1 for line in open(fai_file))
num_contigs_parsed = contig_groups.shape[0] + contig_groups['contigs'][contig_groups.shape[0]-1].count(',')
assert num_contigs_fai == num_contigs_parsed, "Chromosomes in .fai not parsed correctly."
