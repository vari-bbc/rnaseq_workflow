#PBS -l walltime=100:00:00
#PBS -l mem=8gb
#PBS -l nodes=node065:ppn=1
#PBS -m ae
#PBS -N rnaseq_workflow
#PBS -o logs/rnaseq_workflow.o
#PBS -e logs/rnaseq_workflow.e

cd ${PBS_O_WORKDIR}

# export Singularity 3.4.0-1 installed on node065 to $PATH.
# snakemake should now automatically call Singularity 3.4.0-1
PATH=/usr/local/bin/:$PATH

conda activate snakemake
# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --use-conda -n > logs/runs/rnaseq-workflow_${TIME}.txt
snakemake --dag | dot -Tpng > logs/runs/rnaseq-workflow_${TIME}.png

# run snakemake
snakemake \
-s Snakefile \
-j 48 \
--cluster-config src/cluster.yaml \
--cluster 'qsub -q {cluster.qname} -l nodes={cluster.nodes}:ppn={cluster.ppn} -l mem={cluster.mem} -l walltime={cluster.time} -m ea -o error_files/ -e error_files/' \
--use-conda \
--use-singularity

# generate report
#snakemake \
#--report logs/runs/rnaseq-workflow_${TIME}.html
# -o {cluster.std_oe}
