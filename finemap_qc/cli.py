"""Command line interface: finemap-qc --credible-sets cs.tsv.gz [annotations] --out DIR"""
import argparse
import os
import sys

import pandas as pd

from . import annotate, io, report, stats, uncertainty


def _named(value):
    if "=" not in value:
        raise argparse.ArgumentTypeError(f"expected NAME=PATH, got {value!r}")
    name, path = value.split("=", 1)
    return name, path


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="finemap-qc",
        description="Enrichment of high-PIP vs low-PIP variants for features indicative of causality.")
    p.add_argument("--credible-sets", required=True,
                   help="TSV with cs_id, variant (chr_pos_ref_alt, GRCh38), pip")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--gtf", help="GENCODE GTF: adds distance to nearest protein-coding TSS and splice site")
    p.add_argument("--canonical-only", action="store_true",
                   help="use only Ensembl_canonical transcripts from the GTF")
    p.add_argument("--bed", action="append", type=_named, default=[], metavar="NAME=PATH",
                   help="interval feature, e.g. encode_ccre=GRCh38-cCREs.bed (repeatable)")
    p.add_argument("--variant-set", action="append", type=_named, default=[], metavar="NAME=PATH",
                   help="variant ID list feature, e.g. missense=missense.txt (repeatable)")
    p.add_argument("--score", action="append", type=_named, default=[], metavar="NAME=PATH",
                   help="quantitative score table (variant, score), e.g. spliceai=spliceai.tsv (repeatable)")
    p.add_argument("--high-pip", type=float, default=0.8)
    p.add_argument("--low-pip", type=float, default=0.05)
    p.add_argument("--anchor-pip", type=float, default=0.5,
                   help="low_anchored: low-PIP variants from credible sets whose top PIP exceeds this")
    p.add_argument("--duplicate-policy", choices=["error", "collapse-identical"], default="error")
    p.add_argument("--bootstrap-replicates", type=int, default=1000)
    p.add_argument("--bootstrap-seed", type=int, default=0)
    args = p.parse_args(argv)
    if not (0 <= args.low_pip < args.high_pip <= 1 and 0 <= args.anchor_pip <= 1):
        p.error("thresholds must satisfy 0 <= low < high <= 1 and 0 <= anchor <= 1")
    if args.bootstrap_replicates < 100:
        p.error("--bootstrap-replicates must be >= 100")
    return args


def annotate_variants(variants, args):
    """Return (annotation DataFrame indexed by variant, binary feature names, quantitative feature names)."""
    annot = variants.set_index("variant")
    chroms, pos = annot["chrom"].to_numpy(), annot["pos"].to_numpy()
    binary, quant = [], []
    if args.gtf:
        genes = annotate.load_gene_features(args.gtf, args.canonical_only)
        annot["dist_tss"] = annotate.nearest_distance(chroms, pos, genes["tss"])
        annot["dist_splice_site"] = annotate.nearest_distance(chroms, pos, genes["splice"])
        quant += ["dist_tss", "dist_splice_site"]
    for name, path in args.bed:
        annot[name] = annotate.overlaps_intervals(chroms, pos, io.read_bed(path))
        binary.append(name)
    for name, path in args.variant_set:
        annot[name] = annot.index.isin(io.read_variant_set(path))
        binary.append(name)
    for name, path in args.score:
        annot[name] = io.read_scores(path).reindex(annot.index).to_numpy()
        quant.append(name)
    names = binary + quant
    if len(set(names)) != len(names):
        sys.exit("feature names must be unique")
    return annot, binary, quant


def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    out = lambda f: os.path.join(args.out, f)  # noqa: E731

    cs = io.read_credible_sets(args.credible_sets, args.duplicate_policy)
    variants = cs[["variant", "chrom", "pos"]].drop_duplicates("variant")
    annot, binary, quant = annotate_variants(variants, args)

    groups = stats.define_groups(cs, args.high_pip, args.low_pip, args.anchor_pip)
    summary, per_cs = stats.describe(cs, groups, args.anchor_pip)
    summary["duplicate_rows_removed"] = cs.attrs.get("duplicate_rows_removed", 0)
    bin_tab = stats.binary_enrichment(annot, groups, binary)
    q_tab = stats.quantitative_comparison(annot, groups, quant)

    clusters = uncertainty.variant_clusters(cs)
    bin_tab = uncertainty.cluster_intervals(bin_tab, annot, groups, clusters, True,
                                           args.bootstrap_replicates, args.bootstrap_seed)
    q_tab = uncertainty.cluster_intervals(q_tab, annot, groups, clusters, False,
                                         args.bootstrap_replicates, args.bootstrap_seed)
    clusters.to_csv(out("variant_clusters.tsv"), sep="\t", index_label="variant")

    summary.to_frame("value").to_csv(out("descriptive_summary.tsv"), sep="\t", index_label="metric")
    per_cs.to_csv(out("credible_sets.tsv"), sep="\t")
    bin_tab.to_csv(out("binary_enrichment.tsv"), sep="\t", index=False)
    q_tab.to_csv(out("quantitative_comparison.tsv"), sep="\t", index=False)
    group_col = pd.Series("other", index=annot.index)
    for g in ("low_all", "low_anchored", "high"):  # most specific label wins
        group_col[group_col.index.isin(groups[g])] = g
    annot.assign(group=group_col).to_csv(out("variant_annotations.tsv.gz"), sep="\t")

    images = {"descriptive": [out("cs_sizes.png")], "binary": [], "quantitative": []}
    report.plot_cs_sizes(per_cs, images["descriptive"][0])
    if binary:
        images["binary"] += [out("forest.png"), out("binary_features.png")]
        report.plot_forest(bin_tab, images["binary"][0])
        report.plot_binary(bin_tab, images["binary"][1])
    for f in quant:
        path = out(f"hist_{f}.png")
        report.plot_quantitative(annot, groups, f, path, log_distance=f.startswith("dist_"))
        images["quantitative"].append(path)

    name = lambda p: p and os.path.basename(p)  # noqa: E731 - reports get shared; don't leak local paths
    params = {"credible_sets": name(args.credible_sets), "high_pip": f"> {args.high_pip}",
              "low_pip": f"< {args.low_pip}", "anchor_pip (low_anchored)": f"> {args.anchor_pip}",
              "gtf": name(args.gtf), "canonical_only": args.canonical_only,
              "duplicate_policy": args.duplicate_policy,
              "duplicate_rows_removed": cs.attrs.get("duplicate_rows_removed", 0),
              "ci_method": "cluster percentile bootstrap (95%)",
              "cluster_unit": "locus_id + shared variants" if "locus_id" in cs else "credible sets + shared variants",
              "bootstrap_replicates": args.bootstrap_replicates, "bootstrap_seed": args.bootstrap_seed}
    params.update({f"bed:{n}": name(p) for n, p in args.bed})
    params.update({f"variant_set:{n}": name(p) for n, p in args.variant_set})
    params.update({f"score:{n}": name(p) for n, p in args.score})
    report.write_html(out("report.html"), "Fine-mapping QC report", params, summary, bin_tab, q_tab, images)
    print(f"Report written to {out('report.html')}")


if __name__ == "__main__":
    main()
