rule addRG:
    input:
        "analysis/star/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bam="analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam",
    log:
        out="logs/00_addRG/{sample}.o",
        err="logs/00_addRG/{sample}.e"
    benchmark:
        "benchmarks/00_addRG/{sample}.txt"
    envmodules:
        "bbc/picard/picard-2.23.3"
    params:
        unit=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['unit'].values[0],
        lib=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['library'].values[0],
        sample=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['sample'].values[0],
        platform_unit=lambda wildcards: var_calling_units[var_calling_units["unit"]==wildcards.sample]['platform_unit'].values[0]
    threads: 4
    resources: 
        mem_gb = 64
    shell:
        """
        java -Xms16g -Xmx64g -Djava.io.tmpdir=./tmp -jar $PICARD AddOrReplaceReadGroups \
        --INPUT {input} \
        --OUTPUT {output.bam} \
        --RGID {params.unit} \
        --RGLB {params.lib} \
        --RGPL ILLUMINA \
        --RGPU {params.platform_unit} \
        --RGSM {params.sample} 1>{log.out} 2>{log.err}
        """


# variant calling based on the GATK best practices as documented at https://github.com/gatk-workflows/gatk4-rnaseq-germline-snps-indels/blob/master/gatk4-rna-best-practices.wdl -- accessed Aug 11, 2020

rule markdups:
    input:
        lambda wildcards: expand("analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam", sample=var_calling_units[var_calling_units["sample"]==wildcards.sample].index.values)
    output:
        bam="analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.bam",
        metrics="analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.metrics"
    params: 
        input_params=lambda wildcards: expand("--INPUT analysis/00_addRG/{sample}.Aligned.sortedByCoord.out.addRG.bam", sample=var_calling_units[var_calling_units["sample"]==wildcards.sample].index.values)
    log:
        out="logs/01_markdup/{sample}.o",
        err="logs/01_markdup/{sample}.e"
    benchmark:
        "benchmarks/01_markdup/{sample}.txt"
    envmodules:
        "bbc/picard/picard-2.23.3"
    threads: 4
    resources: 
        mem_gb = 64
    shell:
        """
        java -Xms16g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $PICARD MarkDuplicates \
        {params.input_params} \
        --OUTPUT {output.bam} \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY SILENT \
        --METRICS_FILE {output.metrics} 1>{log.out} 2>{log.err}
        """

rule splitncigar:
    input:
        "analysis/01_markdup/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.bam"
    output:
        "analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam"
    log:
        out="logs/02_splitncigar/{sample}.o",
        err="logs/02_splitncigar/{sample}.e"
    benchmark:
        "benchmarks/02_splitncigar/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 64
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SplitNCigarReads \
        -R {params.ref_fasta} \
        -I {input} \
        -O {output} 1>{log.out} 2>{log.err}

        """

rule base_recalibrate:
    input:
        "analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam"
    output:
        "analysis/03_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam.recal_data.table"
    log:
        out="logs/03_base_recal/{sample}.o",
        err="logs/03_base_recal/{sample}.e"
    benchmark:
        "benchmarks/03_base_recal/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        known_snps=config["ref"]["known_snps"],
        known_indels=config["ref"]["known_indels"]
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 64
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=./tmp" \
            BaseRecalibrator \
            -R {params.ref_fasta} \
            -I {input} \
            -O {output} \
            -known-sites {params.known_snps} \
            -known-sites {params.known_indels} 1>{log.out} 2>{log.err}

        """

rule applyBQSR:
    input:
        bam="analysis/02_splitncigar/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam",
        recal_table="analysis/03_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.bam.recal_data.table"
    output:
        "analysis/04_apply_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.bam"
    log:
        out="logs/04_apply_base_recal/{sample}.o",
        err="logs/04_apply_base_recal/{sample}.e"
    benchmark:
        "benchmarks/04_apply_base_recal/{sample}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 64
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g  -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -Djava.io.tmpdir=./tmp" \
            ApplyBQSR \
            --add-output-sam-program-record \
            -R {params.ref_fasta} \
            -I {input.bam} \
            -O {output} \
            --bqsr-recal-file {input.recal_table} 1>{log.out} 2>{log.err}

        """

rule haplotypecaller:
    input:
        bam="analysis/04_apply_base_recal/{sample}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.bam"
    output:
        "analysis/05_haplotypecaller/{sample}.{contig_group}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.g.vcf.gz"
    log:
        out="logs/05_haplotypecaller/{sample}.{contig_group}.o",
        err="logs/05_haplotypecaller/{sample}.{contig_group}.e"
    benchmark:
        "benchmarks/05_haplotypecaller/{sample}.{contig_group}.txt"
    params:
        dbsnp=config["ref"]["known_snps"],
        ref_fasta=config["ref"]["sequence"],
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
        
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 80
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
        --dbsnp {params.dbsnp} \
        {params.contigs} 1>{log.out} 2>{log.err}

        """

rule combinevar:
    input:
        lambda wildcards: expand("analysis/05_haplotypecaller/{sample}.{contig_group}.Aligned.sortedByCoord.out.addRG.mrkdup.splitncigar.baserecal.g.vcf.gz", sample=var_calling_units['sample'].unique(), contig_group=wildcards.contig_group)

    output:
        touch=touch("analysis/06_combinevar/{contig_group}.done"),
        genomicsdb=directory("analysis/06_combinevar/{contig_group}.genomicsdb"),
    log:
        out="logs/06_combinevar/all.{contig_group}.o",
        err="logs/06_combinevar/all.{contig_group}.e"
    benchmark:
        "benchmarks/06_combinevar/{contig_group}.txt"
    params:
        sample_gvcfs = lambda wildcards, input: list(map("-V {}".format, input)),
        contigs = lambda wildcards: "-L " + contig_groups[contig_groups.name == wildcards.contig_group]['contigs'].values[0].replace(",", " -L "),
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenomicsDBImport \
        {params.sample_gvcfs} \
        --genomicsdb-workspace-path {output.genomicsdb} \
        {params.contigs} 1>{log.out} 2>{log.err}

        """

rule jointgeno:
    input:
        "analysis/06_combinevar/{contig_group}.done"
        #"analysis/06_combinevar/{contig_group}.genomicsdb"
    output:
        vcf="analysis/07_jointgeno/all.{contig_group}.vcf.gz",
    log:
        out="logs/07_jointgeno/all.{contig_group}.o",
        err="logs/07_jointgeno/all.{contig_group}.e"
    benchmark:
        "benchmarks/07_jointgeno/all.{contig_group}.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        genomicsdb="analysis/06_combinevar/{contig_group}.genomicsdb"
    envmodules:
        "bbc/gatk/gatk-4.1.8.1"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        GenotypeGVCFs \
        -R {params.ref_fasta} \
        -V gendb://{params.genomicsdb} \
        -O {output.vcf} 1>>{log.out} 2>>{log.err}
        """

rule merge_and_filter_vcf:
    input:
        expand("analysis/07_jointgeno/all.{contig_grp}.vcf.gz", contig_grp=contig_groups.name)
    output:
        raw="analysis/08_merge_and_filter/all.merged.vcf.gz",
        filt="analysis/08_merge_and_filter/all.merged.filt.vcf.gz",
        pass_only="analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz",
        vt_peek_raw="analysis/08_merge_and_filter/all.merged.vcf.gz.vt_peek.txt",
        vt_peek_pass="analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz.vt_peek.txt"
    log:
        out="logs/08_merge_and_filter/out.o",
        err="logs/08_merge_and_filter/out.e"
    benchmark:
        "benchmarks/08_merge_and_filter/benchmark.txt"
    params:
        ref_fasta=config["ref"]["sequence"],
        in_vcfs=expand("--INPUT analysis/07_jointgeno/all.{contig_grp}.vcf.gz", contig_grp=contig_groups.name)
    envmodules:
        "bbc/gatk/gatk-4.1.8.1",
        "bbc/vt/vt-0.1.16"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        MergeVcfs \
        {params.in_vcfs} \
        --OUTPUT {output.raw} \
        1>>{log.out} 2>>{log.err}

        vt peek -r {params.ref_fasta} {output.raw} 2> {output.vt_peek_raw} 1>>{log.out}

        echo "mergeVcfs done." >> {log.out}
        echo "mergeVcfs done." >> {log.err}

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        VariantFiltration \
        --R {params.ref_fasta} \
        --V {output.raw} \
        --window 35 \
        --cluster 3 \
        --filter-name "FS" \
        --filter "FS > 30.0" \
        --filter-name "QD" \
        --filter "QD < 2.0" \
        -O {output.filt} \
        1>>{log.out} 2>>{log.err}
        
        echo "VariantFiltration done." >> {log.out}
        echo "VariantFiltration done." >> {log.err}

        gatk --java-options "-Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp" \
        SelectVariants \
        -R {params.ref_fasta} \
        -V {output.filt} \
        --exclude-filtered \
        -O {output.pass_only} \
        1>>{log.out} 2>>{log.err}
        
        echo "SelectVariants done." >> {log.out}
        echo "SelectVariants done." >> {log.err}

        vt peek -r {params.ref_fasta} {output.pass_only} 2> {output.vt_peek_pass} 1>>{log.out}
        """

rule variant_annot:
    input:
        "analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        html="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.html",
        vcf="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz",
        tbi="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff.vcf.gz.tbi",
        html_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.html",
        vcf_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz",
        tbi_canon="analysis/09a_variant_annot/all.merged.filt.PASS.snpeff_canonical.vcf.gz.tbi",
    log:
        out="logs/09a_variant_annot/out.o",
        err="logs/09a_variant_annot/out.e"
    benchmark:
        "benchmarks/09a_variant_annot/benchmark.txt"
    params:
        db_id=config["ref"]["snpeff_db_id"],
    envmodules:
        "bbc/SnpEff/SnpEff-4.3t",
        "bbc/htslib/htslib-1.10.2"
    threads: 4
    resources: 
        mem_gb = 80
    shell:
        """
        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -canon \
        -onlyProtein \
        -stats {output.html_canon} \
        {params.db_id} \
        {input} \
        2>>{log.err} | \
        bgzip > {output.vcf_canon}

        tabix {output.vcf_canon} 2>>{log.err} 1>>{log.out}


        java -Xms8g -Xmx{resources.mem_gb}g -Djava.io.tmpdir=./tmp -jar $SNPEFF/snpEff.jar eff \
        -v \
        -onlyProtein \
        -stats {output.html} \
        {params.db_id} \
        {input} \
        2>>{log.err} | \
        bgzip > {output.vcf}

        tabix {output.vcf} 2>>{log.err} 1>>{log.out}
        """


rule snprelate:
    input:
        "analysis/08_merge_and_filter/all.merged.filt.PASS.vcf.gz"
    output:
        "analysis/09b_snp_pca_and_dendro/report.html"
    params:
        gds="analysis/09b_snp_pca_and_dendro/all.gds"
    conda:
        "envs/R.yaml"
    #envmodules:
    #    "bbc/cairo/cairo-1.16.0",
    #    "bbc/R/R-3.6.0",
    #    "bbc/pandoc/pandoc-2.7.3",
    threads: 1
    resources:
        mem_gb = 60
    script:
        "bin/snprelate.Rmd"
