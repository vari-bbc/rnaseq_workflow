rule get_mapping_rates:
    input:
        # forces count_matrix to run first
        "deliverables/counts.tsv",
    output:
        rates = "deliverables/UniquelyMappingRates.txt",
        reads = "deliverables/UniquelyMappingReads.txt",
    shell:
        "bash src/alignment_stats.sh {output.rates} {output.reads}"
