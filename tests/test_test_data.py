import csv
import gzip


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_test_config_uses_repo_local_fixtures(repo_root, test_config):
    assert test_config["quick_ref"]["species_name"] is None
    assert test_config["PE_or_SE"] == "PE"
    assert test_config["call_variants"] is False
    assert test_config["run_vis_bigwig"] is False
    assert test_config["run_rseqc"] is False
    assert test_config["deploy_to_shinyio"] is False

    fixture_paths = [
        test_config["ref"]["annotation"],
        test_config["ref"]["dict"],
        test_config["ref"]["fai"],
        test_config["ref"]["known_indels"],
        test_config["ref"]["known_snps"],
        test_config["ref"]["sequence"],
        test_config["ref"]["index"],
        test_config["ref"]["salmon_index"],
        test_config["grouped_contigs"],
        test_config["units"],
        test_config["comparisons"],
    ]
    for fixture in fixture_paths:
        assert (repo_root / fixture).exists(), fixture


def test_test_units_match_fastq_fixtures(repo_root, test_config):
    units = read_tsv(repo_root / test_config["units"])
    assert units
    assert {"sample", "group", "fq1", "fq2", "RG"} <= set(units[0])

    seen_fastqs = set()
    for row in units:
        assert row["sample"]
        assert row["group"]
        for column in ("fq1", "fq2"):
            fastq = row[column]
            fastq_path = repo_root / "raw_data" / fastq
            assert fastq_path.exists(), fastq
            assert fastq not in seen_fastqs
            seen_fastqs.add(fastq)
            with gzip.open(fastq_path, "rt") as handle:
                assert handle.readline().startswith("@")

    group_counts = {}
    for row in units:
        group_counts[row["group"]] = group_counts.get(row["group"], 0) + 1
    assert all(count >= 2 for count in group_counts.values())


def test_test_comparisons_match_units(repo_root, test_config):
    units = read_tsv(repo_root / test_config["units"])
    comparisons = read_tsv(repo_root / test_config["comparisons"])
    groups = {row["group"] for row in units}

    assert comparisons
    assert {"comparison_name", "group_test", "group_reference", "group_reg_formula"} <= set(comparisons[0])
    for row in comparisons:
        assert row["group_test"] in groups
        assert row["group_reference"] in groups
        assert row["group_reg_formula"].startswith("~")


def test_grouped_contigs_match_reference_fai(repo_root, test_config):
    contigs = read_tsv(repo_root / test_config["grouped_contigs"])
    fai_contigs = {
        line.split("\t", 1)[0]
        for line in (repo_root / test_config["ref"]["fai"]).read_text().splitlines()
        if line
    }
    grouped_contigs = set()
    for row in contigs:
        grouped_contigs.update(contig for contig in row["contigs"].split(",") if contig)

    assert grouped_contigs == fai_contigs
