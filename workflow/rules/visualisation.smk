def get_bigwig_norm_factor(wildcards, input):
    df = pd.read_table(input.norm_factors)
    scalefactor = df[df['sample']==wildcards.sample]['sizefactor'].values[0]
    return str(scalefactor)

def get_bw_strand_param (wildcards):
    if wildcards.dir=="fwd":
        return "--filterRNAstrand forward"
    elif wildcards.dir=="rev":
        return "--filterRNAstrand reverse"
    else:
        return " "

rule bigwigs:
    input:
        bam="results/star/{sample}.sorted.bam",
        norm_factors="results/SummarizedExperiment/DESeq2_sizeFactors_reciprocal.tsv"
    output:
        bw="results/bigwigs/{sample}.{dir}.bw"
    params:
        scale_factor=get_bigwig_norm_factor,
        strand=get_bw_strand_param,
    benchmark:
        "benchmarks/bigwigs/{sample}.{dir}.txt"
    envmodules:
        config['modules']['deeptools']
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 16,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        bamCoverage -b {input.bam} -o {output.bw}  --minMappingQuality 30  \
                --normalizeUsing "None" -p {threads} --binSize 10 \
                --scaleFactor {params.scale_factor} {params.strand}
        """

rule avg_bigwigs:
    input:
        lambda wildcards: expand("results/bigwigs/{sample}.{dir}.bw", sample=samples[samples['group']==wildcards.group]['sample'], dir=wildcards.dir)
    output:
        bw="results/avg_bigwigs/{group}.{dir}.bw"
    params:
    benchmark:
        "benchmarks/avg_bigwigs/{group}.{dir}.txt"
    envmodules:
        config['modules']['deeptools']
    threads: 8
    resources:
        nodes = 1,
        mem_gb = 72,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        bigwigAverage -b {input} --binSize 10 -p {threads} -o {output.bw} -of "bigwig"
        """

