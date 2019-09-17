#PBS -l walltime=100:00:00
#PBS -l mem=8gb
#PBS -l nodes=1:ppn=1
#PBS -M dean.pettinga@vai.org
#PBS -m ae
#PBS -N rnaseq_workflow
#PBS -o logs/rnaseq_workflow.o
#PBS -e logs/rnaseq_workflow.o

cd ${PBS_O_WORKDIR}
# save DAG job file with time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")
snakemake --use-conda -n > logs/runs/rnaseq-workflow_${TIME}.txt
snakemake --dag | dot -Tpng > logs/runs/rnaseq-workflow_${TIME}.png

snakemake \
-s Snakefile \
-j 48 \
--cluster-config src/cluster.yaml \
--cluster 'qsub -q {cluster.qname} -l nodes={cluster.nodes}:ppn={cluster.ppn} -l mem={cluster.mem} -l walltime={cluster.time} -m ea -o error_files/ -e error_files/' \
--use-conda \

# generate report
#snakemake \
#--report logs/runs/rnaseq-workflow_${TIME}.html
# -o {cluster.std_oe}
