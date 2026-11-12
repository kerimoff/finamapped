#!/usr/bin/env python3
"""Build a list of 'replicating' fine-mapped variants for use with finemap-qc --variant-set.

A variant is kept if it has PIP > --pip in at least --min-studies explicitly
identified studies. Supply one --study-ids value per file, repeating IDs for
tissues/traits from the same cohort. Different IDs do not prove independence.

Example:
    build_replicating_set.py a.tsv b.tsv --study-ids cohortA cohortB --out rep.txt

FinnGen inputs must use coordinate IDs, not rsIDs. Select the actual coordinate
variant and PIP columns with --variant-col and --pip-col.

Do not include the dataset you are evaluating, or the enrichment will be circular.
"""
import argparse
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from finemap_qc.references import high_pip_variants  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+")
    p.add_argument("--out", required=True)
    p.add_argument("--pip", type=float, default=0.8)
    p.add_argument("--min-studies", type=int, default=2)
    p.add_argument("--variant-col", default="variant")
    p.add_argument("--pip-col", default="pip")
    p.add_argument("--study-ids", nargs="+", required=True,
                   help="one independent study/cohort ID per input file; repeated IDs are pooled")
    p.add_argument("--exclude-study", nargs="*", default=[])
    a = p.parse_args(argv)
    if len(a.study_ids) != len(a.files):
        p.error("provide exactly one --study-ids value per input file")
    if not 0 <= a.pip <= 1 or a.min_studies < 2:
        p.error("require 0 <= --pip <= 1 and --min-studies >= 2")

    counts = Counter()
    study_hits = {}
    for f, study in zip(a.files, a.study_ids):
        if study in a.exclude_study:
            continue
        df = pd.read_csv(f, sep="\t", usecols=[a.variant_col, a.pip_col], dtype={a.variant_col: str},
                         compression="gzip" if f.endswith((".gz", ".bgz")) else "infer")
        study_hits.setdefault(study, set()).update(
            high_pip_variants(df, a.pip, a.variant_col, a.pip_col))
    for hits in study_hits.values():
        counts.update(hits)
    keep = sorted(v for v, n in counts.items() if n >= a.min_studies)
    Path(a.out).write_text("".join(v + "\n" for v in keep))
    print(f"{len(keep)} variants with PIP > {a.pip} in >= {a.min_studies} of {len(study_hits)} studies -> {a.out}")


if __name__ == "__main__":
    main()
