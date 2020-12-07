#PBS -l walltime=100:00:00
#PBS -l mem=8gb
#PBS -l nodes=node065
#PBS -m ae
#PBS -N rnaseq_workflow
#PBS -o logs/rnaseq_workflow.o
#PBS -e logs/rnaseq_workflow.e

cd ${PBS_O_WORKDIR}
# export Singularity 3.4.0-1 installed on node065 to $PATH.
# snakemake should now automatically call Singularity 3.4.0-1
#PATH=/usr/local/bin/:$PATH

module load bbc/snakemake/snakemake-5.28.0

# specify which conda installation to use
conda_setup='/primary/vari/software/BBC/miniconda_bare/etc/profile.d/conda.sh'

# this make 'conda' callable and allows conda environments to be created.
source $conda_setup

# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")

# make logs dir if it does not exist already. Without this, logs/ is automatically generate only after the first run of the pipeline
logs_dir="logs/runs"
[[ -d $logs_dir ]] || mkdir -p $logs_dir

snakemake --use-envmodules -n > logs/runs/workflow_${TIME}.txt
snakemake --dag | dot -Tpng > logs/runs/workflow_${TIME}.png


snakemake \
-p \
--use-envmodules \
--use-conda \
--jobs 100 \
--cluster "qsub \
-q bbc \
-V \
-l nodes=1:ppn={threads} \
-l mem={resources.mem_gb}gb \
-l walltime=100:00:00 \
-o logs/runs/ \
-e logs/runs/"
