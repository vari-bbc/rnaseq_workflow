def get_strandness(samples):
    if "strandedness" in samples.columns:
        return samples["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*samples.shape[0]

rule mapping_summary:
    # get unique mapping rates to text file from STAR logs
    input:
        expand("star/{samples.sample}.ReadsPerGene.out.tab", samples=samples.itertuples())
    output:
        rates = "deliverables/UniquelyMappingRates.txt",
        reads = "deliverables/UniquelyMappingReads.txt",
        matrix = "deliverables/starMatrix.txt"
    log:
        "logs/mapping_summary.log"
    shell:
        "bash src/get_star_matrix.sh {output.rates} {output.reads} {output.matrix}"
