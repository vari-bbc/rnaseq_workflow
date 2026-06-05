import re


RULE_RE = re.compile(r"^(?:rule|checkpoint)\s+([A-Za-z_][A-Za-z0-9_]*)\s*:", re.MULTILINE)
INLINE_INCLUDE_RE = re.compile(r"^\s*include:\s*[\"']([^\"']+)[\"']")
QUOTED_PATH_RE = re.compile(r"^[\"']([^\"']+)[\"']")


def snakefile_includes(text):
    includes = []
    lines = text.splitlines()
    for index, line in enumerate(lines):
        inline_match = INLINE_INCLUDE_RE.match(line)
        if inline_match:
            includes.append(inline_match.group(1))
            continue

        if line.strip() != "include:":
            continue

        for next_line in lines[index + 1:]:
            stripped = next_line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            match = QUOTED_PATH_RE.match(stripped)
            if match:
                includes.append(match.group(1))
            break
    return includes


def test_expected_rules_are_defined(repo_root):
    rule_files = [repo_root / "workflow/Snakefile", *sorted((repo_root / "workflow/rules").glob("*.smk"))]
    rules = set()
    for path in rule_files:
        rules.update(RULE_RE.findall(path.read_text()))

    expected_rules = {
        "all",
        "add_DE_to_SE",
        "avg_bigwigs",
        "bigwigs",
        "BQSR",
        "CollectRnaSeqMetrics",
        "combinevar",
        "concat_fastqs",
        "deploy_isee_to_shinyappio",
        "deseq2",
        "fastq_screen",
        "fastqc",
        "filter_vcf",
        "get_rRNA_intervals_from_gtf",
        "gsea",
        "haplotypecaller",
        "isee",
        "jointgeno",
        "make_final_report",
        "make_genes_bed",
        "make_genes_ref_flat",
        "make_Rproject",
        "markdups",
        "merge_vcf",
        "multiqc",
        "rename_fastqs",
        "rseqc_genebody_cov",
        "salmon",
        "seqtk",
        "snprelate",
        "sortmerna",
        "sortVCF",
        "splitncigar",
        "STAR",
        "SummarizedExperiment",
        "trim_galore_PE",
        "trim_galore_SE",
        "variant_annot",
    }
    assert expected_rules <= rules


def test_snakefile_includes_existing_rule_files(repo_root):
    snakefile = repo_root / "workflow/Snakefile"
    includes = snakefile_includes(snakefile.read_text())
    assert includes
    for include in includes:
        assert (snakefile.parent / include).exists(), include
