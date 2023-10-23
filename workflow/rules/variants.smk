# variant calling based on the GATK best practices as documented at https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl -- accessed Aug 11, 2020

rule markdups:
    input:
        "results/star/{sample}.Aligned.out.bam"
    output:
        bam="results/variant_calling/markdup/{sample}.Aligned.out.mrkdup.bam",
        metrics="results/variant_calling/markdup/{sample}.Aligned.out.mrkdup.metrics"
    params: 
    benchmark:
        "benchmarks/variant_calling/markdup/{sample}.txt"
    envmodules:
        config['modules']['gatk']
    threads: 4
    resources: 
        mem_gb = 120,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
            MarkDuplicatesSpark \
            --spark-master local[{threads}] \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            --conf spark.executor.cores={threads} \
            --conf spark.local.dir=./tmp \
            --conf spark.driver.memory=6g \
            --conf spark.executor.memory=5g

        """

rule splitncigar:
    input:
        "results/variant_calling/markdup/{sample}.Aligned.out.mrkdup.bam"
    output:
        "results/variant_calling/splitncigar/{sample}.Aligned.out.mrkdup.splitncigar.bam"
    benchmark:
        "benchmarks/variant_calling/splitncigar/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        config['modules']['gatk']
    threads: 4
    resources: 
        mem_gb = 96,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SplitNCigarReads \
        -R {params.ref_fasta} \
        -I {input} \
        -O {output} 
        """

rule haplotypecaller:
    input:
        bam=lambda wildcards: "results/variant_calling/splitncigar/{sample}.Aligned.out.mrkdup.splitncigar.bam" if wildcards.round == "bqsr" else "results/variant_calling/final/00_BQSR/{sample}.bqsr.bam"
    output:
        "results/variant_calling/{round}/01_haplotypecaller/{sample}.{contig_group}.g.vcf.gz"
    benchmark:
        "benchmarks/variant_calling/{round}/01_haplotypecaller/{sample}.{contig_group}.txt"
    params:
        dbsnp=lambda wildcards: f'--dbsnp {config["ref"]["known_snps"]}' if config["ref"]["known_snps"] != "" else "",
        ref_fasta=config["ref"]["sequence"],
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
        
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        HaplotypeCaller \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -O {output} \
        -ERC GVCF \
        -dont-use-soft-clipped-bases \
        --native-pair-hmm-threads {threads} \
        --standard-min-confidence-threshold-for-calling 20 \
        {params.dbsnp} \
        {params.contigs}
        """

rule combinevar:
    input:
        lambda wildcards: expand("results/variant_calling/{{round}}/01_haplotypecaller/{sample}.{contig_group}.g.vcf.gz", sample=samples['sample'].unique(), contig_group=wildcards.contig_group)

    output:
        touch=touch("results/variant_calling/{round}/02_combinevar/{contig_group}.done"),
        genomicsdb=directory("results/variant_calling/{round}/02_combinevar/{contig_group}.genomicsdb"),
    benchmark:
        "benchmarks/variant_calling/{round}/02_combinevar/{contig_group}.txt"
    params:
        sample_gvcfs = lambda wildcards, input: list(map("-V {}".format, input)),
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 120,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        {params.sample_gvcfs} \
        --reader-threads {threads} \
        --genomicsdb-workspace-path {output.genomicsdb} \
        {params.contigs}
        """

rule jointgeno:
    input:
        "results/variant_calling/{round}/02_combinevar/{contig_group}.done"
    output:
        vcf="results/variant_calling/{round}/03_jointgeno/all.{contig_group}.vcf.gz",
    benchmark:
        "benchmarks/variant_calling/{round}/03_jointgeno/all.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        genomicsdb="results/variant_calling/{round}/02_combinevar/{contig_group}.genomicsdb"
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {params.ref_fasta} \
        -V gendb://{params.genomicsdb} \
        -O {output.vcf}
        """

rule sortVCF:
    """
    Sort the output VCFs from joint genotyping. Merging errors out sometimes if we do not do this step.
    """
    input:
        vcf="results/variant_calling/{round}/03_jointgeno/all.{contig_group}.vcf.gz",
    output:
        sorted_vcf="results/variant_calling/{round}/04_sortvcf/all.{contig_group}.sort.vcf.gz"
    benchmark:
        "benchmarks/variant_calling/{round}/04_sortvcf/all.{contig_group}.txt"
    params:
        dictionary=config['ref']['dict'],
    envmodules:
        config["modules"]["gatk"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SortVcf \
        -I {input.vcf} \
        -O {output.sorted_vcf} \
        -SD {params.dictionary} 
        """

rule merge_vcf:
    """
    Merge the contig group VCFs into one unified VCF.
    """
    input:
        expand("results/variant_calling/{{round}}/04_sortvcf/all.{contig_grp}.sort.vcf.gz", contig_grp=contig_groups.name)
    output:
        raw="results/variant_calling/{round}/05_merge_vcf/all.merged.vcf.gz",
        raw_tbi="results/variant_calling/{round}/05_merge_vcf/all.merged.vcf.gz.tbi",
        vt_peek_raw="results/variant_calling/{round}/05_merge_vcf/all.merged.vcf.gz.vt_peek.txt",
    benchmark:
        "benchmarks/variant_calling/{round}/05_merge_vcf/benchmark.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        dictionary=config['ref']['dict'],
        in_vcfs = lambda wildcards, input: ' '.join(['--INPUT ' + vcf for vcf in input]) 
    envmodules:
        config["modules"]["gatk"],
        config["modules"]["vt"],
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        MergeVcfs \
        {params.in_vcfs} \
        --SEQUENCE_DICTIONARY {params.dictionary} \
        --OUTPUT {output.raw} 
        
        vt peek -r {params.ref_fasta} {output.raw} 2> {output.vt_peek_raw} 
        """

def get_filt_params (wildcards):
    # parameters adapted from https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
    if (wildcards.var_type == "SNP"):
        return """--cluster-window-size 35 \
        --cluster-size 3 \
        --filter-name 'FS' \
        --filter 'FS > 30.0' \
        --filter-name 'QD' \
        --filter 'QD < 2.0'"""
    elif (wildcards.var_type == "INDEL"):
        return """--cluster-window-size 35 \
        --cluster-size 3 \
        --filter-name 'FS' \
        --filter 'FS > 30.0' \
        --filter-name 'QD' \
        --filter 'QD < 2.0'"""
    else:
        raise Exception("var_type wildcard must be either SNP or INDEL.")

rule filter_vcf:
    """
    Do quality filters. Use different paramters depending on SNPs verus indels ('SNP' or 'INDEL').
    """
    input:
        "results/variant_calling/{round}/05_merge_vcf/all.merged.vcf.gz"
    output:
        raw="results/variant_calling/{round}/06_filter_vcf/all.merged.{var_type}.vcf.gz",
        filt="results/variant_calling/{round}/06_filter_vcf/all.merged.filt.{var_type}.vcf.gz",
        pass_only="results/variant_calling/{round}/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz",
        vt_peek_pass="results/variant_calling/{round}/06_filter_vcf/all.merged.filt.PASS.{var_type}.vcf.gz.vt_peek.txt"
    benchmark:
        "benchmarks/variant_calling/{round}/06_filter_vcf/{var_type}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        filt_params=get_filt_params
    envmodules:
        config["modules"]["gatk"],
        config["modules"]["vt"]
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {input} \
        --select-type-to-include {wildcards.var_type} \
        -O {output.raw}

        echo "SelectVariants 1 done." 
        echo "SelectVariants 1 done." 1>&2
        
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        VariantFiltration \
        --R {params.ref_fasta} \
        --V {output.raw} \
        {params.filt_params} \
        -O {output.filt} 
        
        echo "VariantFiltration done." 
        echo "VariantFiltration done." 1>&2

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {output.filt} \
        --exclude-filtered \
        -O {output.pass_only} 
        
        echo "SelectVariants 2 done." 
        echo "SelectVariants 2 done." 1>&2

        vt peek -r {params.ref_fasta} {output.pass_only} 2> {output.vt_peek_pass} 
        """

rule BQSR:
    """
    Base quality score recalibrator
    """
    input:
        bam="results/variant_calling/splitncigar/{sample}.Aligned.out.mrkdup.splitncigar.bam",
        known_snps_vcf=config["ref"]["known_snps"] if config["ref"]["known_snps"] else "results/variant_calling/bqsr/06_filter_vcf/all.merged.filt.PASS.SNP.vcf.gz",
        known_indels_vcf=config["ref"]["known_indels"] if config["ref"]["known_indels"] else "results/variant_calling/bqsr/06_filter_vcf/all.merged.filt.PASS.INDEL.vcf.gz",
    output:
        recal_table="results/variant_calling/final/00_BQSR/{sample}.bqsr.table",
        recal_bam="results/variant_calling/final/00_BQSR/{sample}.bqsr.bam"
    benchmark:
        "benchmarks/variant_calling/final/00_BQSR/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        config["modules"]["gatk"],
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        BaseRecalibrator \
        -R {params.ref_fasta} \
        -I {input.bam} \
        --known-sites {input.known_snps_vcf} \
        --known-sites {input.known_indels_vcf} \
        -O {output.recal_table} 
 
        echo "BaseRecalibrator done." 
        echo "BaseRecalibrator done." 1>&2
        
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        ApplyBQSR \
        -R {params.ref_fasta} \
        -I {input.bam} \
        -bqsr {output.recal_table} \
        -O {output.recal_bam}

        echo "ApplyBQSR done." 
        echo "ApplyBQSR done." 1>&2

        """

rule variant_annot:
    input:
        "results/variant_calling/final/06_filter_vcf/all.merged.filt.PASS.SNP.vcf.gz"
    output:
        html="results/variant_calling/final/07a_variant_annot/all.merged.filt.PASS.snpeff.html",
        vcf="results/variant_calling/final/07a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz",
        tbi="results/variant_calling/final/07a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi",
        html_canon="results/variant_calling/final/07a_variant_annot/all.merged.filt.PASS.snpeff_canonical.html",
        vcf_canon="results/variant_calling/final/07a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz",
        tbi_canon="results/variant_calling/final/07a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz.tbi",
    benchmark:
        "benchmarks/variant_calling/final/07a_variant_annot/benchmark.txt"
    params:
        db_id=config["ref"]["snpeff_db_id"],
    envmodules:
        config['modules']['snpeff'],
        config['modules']['htslib']
    threads: 4
    resources: 
        mem_gb = 80,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -canon \
        -onlyProtein \
        -stats {output.html_canon} \
        {params.db_id} \
        {input} | \
        bgzip > {output.vcf_canon}

        tabix {output.vcf_canon}


        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -onlyProtein \
        -stats {output.html} \
        {params.db_id} \
        {input} | \
        bgzip > {output.vcf}

        tabix {output.vcf} 
        """

rule snprelate:
    input:
        vcf="results/variant_calling/final/06_filter_vcf/all.merged.filt.PASS.SNP.vcf.gz"
    output:
        symlink_rmd = "results/variant_calling/final/07b_snp_pca_and_dendro/snprelate.Rmd",
        symlink_vcf = "results/variant_calling/final/07b_snp_pca_and_dendro/all.merged.filt.PASS.SNP.vcf.gz",
        html = "results/variant_calling/final/07b_snp_pca_and_dendro/snprelate.html",
        outdir = directory("results/variant_calling/final/07b_snp_pca_and_dendro/snprelate_out_files")
    benchmark:
        "benchmarks/variant_calling/final/07b_snprelate/benchmark.txt"
    params:
        rmd='workflow/scripts/snprelate.Rmd',
        wd = lambda wildcards, output: os.path.dirname(output.html),
        in_vcf = lambda wildcards, input: os.path.basename(input.vcf),
        outdir = lambda wildcards, output: os.path.basename(output.outdir)
    envmodules:
        config['modules']['R']
    threads: 1
    resources:
        mem_gb = 60,
        log_prefix=lambda wildcards: "_".join(wildcards) if len(wildcards) > 0 else "log"
    shell:
        """
        ln -sr {input.vcf} {output.symlink_vcf}
        ln -sr {params.rmd} {output.symlink_rmd}

        cd {params.wd}

        Rscript --vanilla -e "rmarkdown::render('snprelate.Rmd', params = list(in_vcf = '{params.in_vcf}', outdir = '{params.outdir}'))"
        """
