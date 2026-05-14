#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH -J rnaseq_workflow_core_hpc
#SBATCH -o logs/run_core_pipeline_hpc.%j.o
#SBATCH -e logs/run_core_pipeline_hpc.%j.e
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 8
#SBATCH --time 120:00:00
#SBATCH --mem=8G
#SBATCH --partition=long
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

snakemake_module="bbc2/snakemake/snakemake-9.13.2"
if command -v module >/dev/null 2>&1; then
    module load "$snakemake_module"
fi

cores="${SNAKEMAKE_CORES:-8}"
config="${SNAKEMAKE_TEST_CONFIG:-tests/test_config/config_core_hpc.yaml}"

common_args=(
    --cores "${cores}"
    --printshellcmds
    --rerun-incomplete
    --latency-wait 30
)

if [[ -n "${SNAKEMAKE_DEPLOYMENT_ARGS:-}" ]]; then
    # shellcheck disable=SC2206
    deployment_args=(${SNAKEMAKE_DEPLOYMENT_ARGS})
else
    deployment_args=(--use-envmodules)
fi

snakemake \
    "${common_args[@]}" \
    "${deployment_args[@]}" \
    --snakefile workflow/Snakefile \
    --configfile "${config}" \
    results/rename_fastqs/SRR1039508_0_R1.fastq.gz \
    results/trimmed_data/SRR1039508_0_R1_val_1.fq.gz \
    results/trimmed_data/SRR1039508_0_R2_val_2.fq.gz \
    results/concat_fastqs/SRR1039508_R1.fastq.gz \
    results/star/SRR1039508.sorted.bam \
    results/star/SRR1039508.Log.final.out \
    results/salmon/SRR1039508/quant.sf
