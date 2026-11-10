import gzip

import numpy as np
import pandas as pd
import pytest

from finemap_qc import annotate, cli, io, stats

GTF = """\
##comment
chr1\tHAVANA\tgene\t1000\t5000\t.\t+\t.\tgene_id "G1"; gene_type "protein_coding";
chr1\tHAVANA\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "protein_coding"; transcript_type "protein_coding"; tag "Ensembl_canonical";
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t3000\t3100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t4800\t5000\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\ttranscript\t20000\t30000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t20000\t21000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t29000\t30000\t.\t-\t.\tgene_id "G2"; transcript_id "T2"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\ttranscript\t50000\t60000\t.\t+\t.\tgene_id "G3"; transcript_id "T3"; gene_type "lncRNA"; transcript_type "lncRNA";
"""


def test_normalize_variant_ids():
    df = io.normalize_variant_ids(["1_100_A_G", "chr1:100:A:G", "X-5-C-T", "garbage"])
    assert df["variant"].tolist()[:3] == ["chr1_100_A_G", "chr1_100_A_G", "chrX_5_C_T"]
    assert pd.isna(df["variant"].iloc[3])


def test_gene_features(tmp_path):
    p = tmp_path / "a.gtf.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(GTF)
    g = annotate.load_gene_features(p)
    assert g["tss"]["1"].tolist() == [1000, 30000]  # lncRNA excluded, minus strand uses end
    assert g["splice"]["1"].tolist() == [1200, 3000, 3100, 4800, 21000, 29000]
    d = annotate.nearest_distance(np.array(["1", "1", "2"]), np.array([1010, 25000, 5]), g["tss"])
    assert d[0] == 10 and d[1] == 5000 and np.isnan(d[2])
    assert annotate.load_gene_features(p, canonical_only=True)["tss"]["1"].tolist() == [1000]


def test_overlaps_intervals(tmp_path):
    p = tmp_path / "a.bed"
    p.write_text("chr1\t0\t1000\nchr1\t100\t200\nchr1\t5000\t5001\n")
    bed = io.read_bed(p)
    hits = annotate.overlaps_intervals(np.array(["1"] * 5 + ["2"]),
                                       np.array([1, 1000, 1001, 5001, 5002, 10]), bed)
    assert hits.tolist() == [True, True, False, True, False, False]


def test_groups_and_enrichment():
    cs = pd.DataFrame({
        "cs_id": ["A_L1"] * 3 + ["B_L1"] * 3,
        "variant": ["v1", "v2", "v3", "v4", "v5", "v1"],
        "pip": [0.9, 0.02, 0.01, 0.3, 0.03, 0.01],
        "trait": ["A"] * 3 + ["B"] * 3,
    })
    g = stats.define_groups(cs)
    assert g["high"] == {"v1"}
    assert g["low_all"] == {"v2", "v3", "v5"}  # v1 removed: it is high-PIP elsewhere
    assert g["low_anchored"] == {"v2", "v3"}
    summary, per_cs = stats.describe(cs, g)
    assert summary["n_credible_sets"] == 2 and summary["frac_cs_with_pip_gt_0.5"] == 0.5
    annot = pd.DataFrame({"f": [True, True, False, False, False]}, index=["v1", "v2", "v3", "v4", "v5"])
    t = stats.binary_enrichment(annot, g, ["f"]).set_index("comparison")
    assert t.loc["high_vs_low_all", "frac_low"] == pytest.approx(1 / 3)
    assert t.loc["high_vs_low_anchored", "fold_enrichment"] == pytest.approx(2.0)


def test_cli_end_to_end(tmp_path):
    rng = np.random.default_rng(0)
    rows = []
    for i in range(200):
        n = rng.integers(1, 30)
        pips = rng.dirichlet(np.ones(n) * 0.3)
        for j, pip in enumerate(pips):
            pos = 1000 + i * 5000 + j * 50
            rows.append((f"ENSG{i}_L1", f"chr1_{pos}_A_G", pip))
    cs = pd.DataFrame(rows, columns=["cs_id", "variant", "pip"])
    cs.to_csv(tmp_path / "cs.tsv", sep="\t", index=False)
    (tmp_path / "a.gtf").write_text(GTF)
    (tmp_path / "open.bed").write_text("chr1\t0\t500000\n")
    (tmp_path / "set.txt").write_text("variant\n" + "\n".join(cs["variant"].iloc[::3]) + "\n")
    (tmp_path / "score.tsv").write_text("variant\tscore\n" + "\n".join(f"{v}\t{k}" for k, v in enumerate(cs["variant"])))
    out = tmp_path / "out"
    cli.main(["--credible-sets", str(tmp_path / "cs.tsv"), "--out", str(out), "--gtf", str(tmp_path / "a.gtf"),
              "--bed", f"open={tmp_path / 'open.bed'}", "--variant-set", f"s={tmp_path / 'set.txt'}",
              "--score", f"sc={tmp_path / 'score.tsv'}"])
    for f in ["report.html", "binary_enrichment.tsv", "quantitative_comparison.tsv", "hist_dist_tss.png",
              "descriptive_summary.tsv", "variant_annotations.tsv.gz"]:
        assert (out / f).exists(), f
    b = pd.read_csv(out / "binary_enrichment.tsv", sep="\t")
    assert set(b["feature"]) == {"open", "s"} and len(b) == 4
    q = pd.read_csv(out / "quantitative_comparison.tsv", sep="\t")
    assert set(q["feature"]) == {"dist_tss", "dist_splice_site", "sc"}
