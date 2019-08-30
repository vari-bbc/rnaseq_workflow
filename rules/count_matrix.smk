def get_strandness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*units.shape[0]

rule count_matrix:
    input:
        # ask for sorted output to ensure that sort_index_bam occurs first.
        expand("analysis/star/{unit.sample}_{unit.unit}.sorted.bam", unit=units.itertuples()),
        # ReadsPerGene s
        #expand("analysis/star/{unit.sample}_{unit.unit}.ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "deliverables/counts.tsv",
        rates = "deliverables/UniquelyMappingRates.txt",
        reads = "deliverables/UniquelyMappingReads.txt",
        matrix = "deliverables/starMatrix.txt",
    params:
        samples=units["sample"].tolist(),
        strand=get_strandness(units)
    conda:
        "../envs/pandas.yaml"
    log:
        "logs/mapping_summary.log"
    # script:
    #     "../src/count-matrix.py"
    shell:
        """
        bash src/get_star_matrix.sh {output.rates} {output.reads} {output.matrix}
        python ../src/count-matrix.py"
        """
