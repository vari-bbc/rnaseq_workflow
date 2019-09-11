rule get_mapping_rates:
    input:
        # forces count_matrix to run first
        "deliverables/counts.tsv",
    output:
        rates = "deliverables/UniquelyMappingRates.txt",
        reads = "deliverables/UniquelyMappingReads.txt",
        matrix = "deliverables/starMatrix.txt",
    shell:
        "bash src/get_star_matrix.sh {output.rates} {output.reads} {output.matrix}"
