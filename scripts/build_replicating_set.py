#!/usr/bin/env python3
"""Build a list of 'replicating' fine-mapped variants for use with finemap-qc --variant-set.

A variant is kept if it has PIP > --pip in at least --min-studies of the input
files. Each input file counts as one study, so concatenate files that should
count as a single study (e.g. several tissues from the same cohort) beforehand.

Example (eQTL Catalogue gene expression credible sets):
    build_replicating_set.py --pip 0.8 --min-studies 2 --out eqtlcat_ge_rep.txt QTD*.credible_sets.tsv.gz
Example (FinnGen SuSiE *.SUSIE.snp.bgz files):
    build_replicating_set.py --variant-col rsid --pip-col prob --out finngen_rep.txt *.snp.bgz

Do not include the dataset you are evaluating, or the enrichment will be circular.
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from finemap_qc.io import normalize_variant_ids  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+")
    p.add_argument("--out", required=True)
    p.add_argument("--pip", type=float, default=0.8)
    p.add_argument("--min-studies", type=int, default=2)
    p.add_argument("--variant-col", default="variant")
    p.add_argument("--pip-col", default="pip")
    a = p.parse_args()

    counts = Counter()
    for f in a.files:
        df = pd.read_csv(f, sep="\t", usecols=[a.variant_col, a.pip_col], dtype={a.variant_col: str},
                         compression="gzip" if f.endswith((".gz", ".bgz")) else "infer")
        hits = df.loc[pd.to_numeric(df[a.pip_col], errors="coerce") > a.pip, a.variant_col]
        counts.update(set(normalize_variant_ids(hits)["variant"].dropna()))
    keep = sorted(v for v, n in counts.items() if n >= a.min_studies)
    Path(a.out).write_text("".join(v + "\n" for v in keep))
    print(f"{len(keep)} variants with PIP > {a.pip} in >= {a.min_studies} of {len(a.files)} files -> {a.out}")


if __name__ == "__main__":
    main()
