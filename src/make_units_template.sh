#!/bin/bash

set -e
set -u
set -o pipefail

# declare variables for options
raw_data_dir="../raw_reads/"
sample_rgx='^[^_]+'
group_rgx="$sample_rgx"
genotype_rgx="$sample_rgx"
condition_rgx="$sample_rgx"
unit_rgx="$sample_rgx"
strandedness='reverse'
out_file='units_template.tsv'

# help message
print_usage() {
  usage_msg="Usage: ./make_units_template.sh [options]\n"
  usage_msg="${usage_msg}Regex flavor is Perl. If multiple matches, only the first match will be used.\n"
  usage_msg="${usage_msg}PE reads are automatically inferred by having two fastq files matching the sample regex.\n"
  usage_msg="${usage_msg}\t-d raw data directory [${raw_data_dir}]\n"
  usage_msg="${usage_msg}\t-s sample regex [${sample_rgx}]\n"
  usage_msg="${usage_msg}\t-p group regex [${sample_rgx}]\n"
  usage_msg="${usage_msg}\t-g genotype regex [${sample_rgx}]\n"
  usage_msg="${usage_msg}\t-c condition regex [${sample_rgx}]\n"
  usage_msg="${usage_msg}\t-u unit regex [${sample_rgx}]\n"
  usage_msg="${usage_msg}\t-o output file name [${out_file}]\n"
  usage_msg="${usage_msg}\t-x strandedness [${strandedness}]\n"
  usage_msg="${usage_msg}\t-h help page\n"

  printf "$usage_msg"
}

# read in options
while getopts 'd:s:p:g:c:u:o:x:h' flag; do
  case "${flag}" in
    d) raw_data_dir="${OPTARG}" ;; 
    s) sample_rgx="${OPTARG}" ;;
    p) group_rgx="${OPTARG}" ;;
    g) genotype_rgx="${OPTARG}" ;;
    c) condition_rgx="${OPTARG}" ;;
    u) unit_rgx="${OPTARG}" ;;
    o) out_file="${OPTARG}" ;;
    x) strandedness="${OPTARG}" ;;
    h) print_usage
       exit 1 ;;
  esac
done

# get the names of all the fastq files
fastq_files=$(ls "${raw_data_dir}" | grep -P '.*\bfq\b.*|.*\bfastq\b.*') # the \b accounts for for example .fastqc files

# check that output file does not exist
if [ -f ${out_file} ]; then
    echo "Output file ${out_file} exists already."
    exit
fi

# store header in variable to  output to new file at the end
out_file_contents="sample\tgroup\tgenotype\tcondition\tunit\tfq1\tfq2\tstrandedness\n"

# get sample ids based on sample regex
sample_ids=$(echo "$fastq_files" | grep -Po "${sample_rgx}" | sort | uniq)

# loop through each sample
for sample_id in ${sample_ids}
do
    sample_fastqs=$(echo "$fastq_files" | grep -P "$sample_id")
    num_sample_fastqs=$(echo "$sample_fastqs" | wc -l)
    fq1=$(echo "$sample_fastqs" | head -n 1)
    fq2=""

    # check that there are no more than 2 matching files per sample
    if [ "$num_sample_fastqs" -gt 2 ]
    then
        echo "${sample_id} has ${num_sample_fastqs} matching fastq files."
        exit
    elif [ "$num_sample_fastqs" -eq 2 ]
    then
        fq2=$(echo "$sample_fastqs" | tail -n 1)
    fi

    # for each of these, take the first match if there are multiple matches
    sample_group=$(echo "$sample_id" | grep -Po "$group_rgx" | head -n 1)
    sample_genotype=$(echo "$sample_id" | grep -Po "$genotype_rgx" | head -n 1)
    sample_condition=$(echo "$sample_id" | grep -Po "$condition_rgx" | head -n 1)
    sample_unit=$(echo "$sample_id" | grep -Po "$unit_rgx" | head -n 1)

    # assemble the info for this sample
    out_line="${sample_id}\t${sample_group}\t${sample_genotype}\t${sample_condition}\t"
    out_line="${out_line}${sample_unit}\t${fq1}\t${fq2}\t${strandedness}\n"

    # append info for this sample to the output
    out_file_contents="${out_file_contents}${out_line}"
done

# output to output file
printf "$out_file_contents" > $out_file


