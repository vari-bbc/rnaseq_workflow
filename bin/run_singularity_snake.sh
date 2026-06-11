#!/bin/bash
#SBATCH --export=NONE
#SBATCH -J rnaseq_workflow
#SBATCH -o rnaseq_workflow.o
#SBATCH -e rnaseq_workflow.e
#SBATCH --ntasks 1
#SBATCH --time 120:00:00
#SBATCH --mem=8G
#SBATCH --partition=long

cd $SLURM_SUBMIT_DIR

snakemake_module="bbc2/snakemake/snakemake-9.13.2"

module load $snakemake_module

if [[ -n "${SINGULARITY_MODULE:-}" ]]; then
    module load "$SINGULARITY_MODULE"
fi

# make logs dir if it does not exist already.
logs_dir="logs/"
[[ -d $logs_dir ]] || mkdir -p $logs_dir


echo "Start snakemake workflow." >&1
echo "Start snakemake workflow." >&2

singularity_args="${SINGULARITY_ARGS:---bind /varidata:/varidata}"

snakemake \
-p \
--latency-wait 20 \
--use-singularity \
--singularity-args "$singularity_args" \
--jobs 100 \
--executor cluster-generic --cluster-generic-submit-cmd "mkdir -p logs/{rule}; sbatch \
-p ${SLURM_JOB_PARTITION} \
--export=ALL \
--nodes 1 \
--ntasks-per-node {threads} \
--mem={resources.mem_gb}G \
-t 120:00:00 \
-o logs/{rule}/{resources.log_prefix}.o \
-e logs/{rule}/{resources.log_prefix}.e" # SLURM hangs if output dir does not exist, so we create it before running sbatch on the snakemake jobs.

echo "snakemake workflow done." >&1
echo "snakemake workflow done." >&2
