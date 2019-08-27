def get_strandness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*units.shape[0]

rule mapping_summary:
    # get unique mapping rates to text file from STAR logs
    input:
        expand("analysis/star/{units.sample}_{units.unit}.ReadsPerGene.out.tab", units=units.itertuples())
    output:
        rates = "deliverables/UniquelyMappingRates.txt",
        reads = "deliverables/UniquelyMappingReads.txt",
        matrix = "deliverables/starMatrix.txt"
    log:
        "logs/mapping_summary.log"
    shell:
        "bash src/get_star_matrix.sh {output.rates} {output.reads} {output.matrix}"

rule count_matrix:
    input:
        expand("analysis/star/{unit.sample}_{unit.unit}.ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "deliverables/counts.tsv"
    params:
        samples=units["sample"].tolist(),
        strand=get_strandness(units)
    conda:
        "../envs/pandas.yaml"
    script:
        "../src/count-matrix.py"
