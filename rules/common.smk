def is_single_end(sample):
    return pd.isnull(samples.loc[sample, "fq2"])

def get_fastqc_output():
    html = []
    zip = []
    for x in (samples['sample']):
        # if paired data
        if pd.isnull(samples.loc[x,'fq2']):
            html.append(expand("qc/fastqc/{sample}_fastqc.html", sample=x))
            zip.append(expand("qc/fastqc/{sample}_fastqc.zip", sample=x))
        if not pd.isnull(samples.loc[x,'fq2']):
            html.append(expand("qc/fastqc/{sample}_{read}_fastqc.html", read=[1, 2], sample=x))
            zip.append(expand("qc/fastqc/{sample}_{read}_fastqc.zip", read=[1, 2], sample=x))

    flat = lambda l: [item for sublist in l for item in sublist]
    return flat([flat(html),flat(zip)])

def get_fastqc_input():
    reads = []
    for x in (samples['sample']):
        # if paired data
        if pd.isnull(samples.loc[x,'fq2']):
            reads.append(expand("trimmed/{sample}_fastq.gz", sample=x))
        else:
            reads.append(expand("trimmed/{sample}_{read}_fastq.gz", read=[1,2], sample=x))

    flat = lambda l: [item for sublist in l for item in sublist]
    return flat(reads)


# for row in samples.loc["sample"]:
