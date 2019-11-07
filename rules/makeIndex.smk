rule downloadGenomeAndGTF:
    output:
        fasta="ref/gencode/GRCh38.primary_assembly.genome.fa",
        GTF="ref/gencode/gencode.v32.annotation.gtf",
        TE_GTF="ref/GRCh38_rmsk_TE.gtf"
    params:
        fasta="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.primary_assembly.genome.fa.gz",
        GTF="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz",
        TE_GTF="http://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_rmsk_TE.gtf.gz",
    shell:
        """
        wget -O - ${params.fasta} | gunzip -c > `basename ${{{params.fasta}%\.*gz}}`
        wget -O - ${params.GTF} | gunzip -c > `basename ${{{params.GTF}%\.*gz}}`
        wget -O - ${params.TE_GTF} | gunzip -c > `basename ${{{params.TE_GTF}%\.*gz}}`
        """

rule makeIndex:
    input:
        fasta="ref/gencode/GRCh38.primary_assembly.genome.fa",
        GTF="ref/gencode/gencode.v32.annotation.gtf",
    output:
        "ref/gencode/Genome",
        "ref/gencode/Log.out",
        "ref/gencode/SA",
        "ref/gencode/SAindex",
        "ref/gencode/chrLength.txt",
        "ref/gencode/chrName.txt",
        "ref/gencode/chrNameLength.txt",
        "ref/gencode/chrStart.txt",
        "ref/gencode/exonGeTrInfo.tab",
        "ref/gencode/exonInfo.tab",
        "ref/gencode/geneInfo.tab",
        "ref/gencode/genomeParameters.txt",
        "ref/gencode/sjdbInfo.txt",
        "ref/gencode/sjdbList.fromGTF.out.tab",
        "ref/gencode/sjdbList.out.tab",
        "ref/gencode/transcriptInfo.tab",
    params:
        sjdbOverhang=config["readLength"],
        indexDir="ref/gencode/"
    log:
        "logs/makeIndex.log"
    threads: 32
    conda:
        "../envs/STAR.yaml"
    shell:
        """
        STAR \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.GTF} \
        --sjdbOverhang {params.sjdbOverhang} \
        --runMode genomeGenerate \
        --genomeDir {params.indexDir} \
        --runThreadN {threads} \
        2> {log}
        """

rule TEtranscripts_GTF_rename:
    input:
        "ref/GRCh38_rmsk_TE.gtf"
    output:
        "ref/gencode/GRCh38_rmsk_TE_edit.gtf"
    conda:
        "../envs/R.yaml"
    script:
        "TEtranscripts_GTF_rename.Rmd"
