"""Readers for credible sets and annotation files.

All variant IDs are normalised to ``chr{chrom}_{pos}_{ref}_{alt}`` so that
files using ``1_123_A_G``, ``chr1:123:A:G`` or ``1-123-A-G`` all match.
"""
from functools import lru_cache

import numpy as np
import pandas as pd


def normalize_variant_ids(ids):
    """Return a DataFrame with columns variant, chrom, pos (NaN rows for unparsable IDs)."""
    parts = pd.Series(ids, dtype=str).str.split(r"[_:\-]", n=3, regex=True, expand=True)
    parts = parts.reindex(columns=range(4))
    chrom = parts[0].str.replace(r"^chr", "", regex=True)
    pos = pd.to_numeric(parts[1], errors="coerce")
    ok = parts.notna().all(axis=1) & pos.notna()
    variant = ("chr" + chrom + "_" + parts[1] + "_" + parts[2] + "_" + parts[3]).where(ok)
    return pd.DataFrame({"variant": variant, "chrom": chrom.where(ok), "pos": pos.where(ok)})


def _is_number(x):
    try:
        float(x)
        return True
    except (TypeError, ValueError):
        return False


def read_credible_sets(path):
    """Read credible sets.

    Uses columns named cs_id, variant, pip if present (so eQTL Catalogue
    *.credible_sets.tsv.gz files work directly); otherwise the first three
    columns, with or without a header line.  If a molecular_trait_id column
    exists it is used to group credible sets by trait; otherwise the trait is
    the cs_id with a trailing ``_L<number>`` removed (eQTL Catalogue convention).
    """
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str)
    first = raw.iloc[0].tolist()
    if not _is_number(first[2]):  # header line present
        raw = raw.iloc[1:].reset_index(drop=True)
        raw.columns = first
    if {"cs_id", "variant", "pip"} <= set(raw.columns):
        keep = ["cs_id", "variant", "pip"] + (["molecular_trait_id"] if "molecular_trait_id" in raw.columns else [])
        df = raw[keep].copy()
    else:
        df = raw.iloc[:, :3].copy()
        df.columns = ["cs_id", "variant", "pip"]

    df["pip"] = pd.to_numeric(df["pip"], errors="raise")
    parsed = normalize_variant_ids(df["variant"])
    bad = parsed["variant"].isna()
    if bad.any():
        raise ValueError(f"{bad.sum()} unparsable variant IDs, e.g. {df.loc[bad, 'variant'].iloc[0]!r}")
    df["variant"] = parsed["variant"].values
    df["chrom"] = parsed["chrom"].values
    df["pos"] = parsed["pos"].astype("int64").values
    if "molecular_trait_id" in df.columns:
        df = df.rename(columns={"molecular_trait_id": "trait"})
    else:
        df["trait"] = df["cs_id"].str.replace(r"_L\d+$", "", regex=True)
    return df


@lru_cache(maxsize=16)
def read_bed(path):
    """Read the first three columns of a BED file into {chrom: (sorted starts, cumulative max ends)}."""
    bed = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2], comment="#",
                      names=["chrom", "start", "end"], dtype={"chrom": str})
    bed = bed[~bed["chrom"].str.startswith(("track", "browser"))]
    bed["chrom"] = bed["chrom"].str.replace(r"^chr", "", regex=True)
    bed["start"] = bed["start"].astype("int64")
    bed["end"] = bed["end"].astype("int64")
    out = {}
    for chrom, g in bed.sort_values(["chrom", "start"]).groupby("chrom"):
        out[chrom] = (g["start"].to_numpy(), np.maximum.accumulate(g["end"].to_numpy()))
    return out


def read_variant_set(path):
    """Read variant IDs from the first column of a file (header lines are skipped automatically)."""
    ids = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0], dtype=str, comment="#")[0]
    return set(normalize_variant_ids(ids)["variant"].dropna())


def read_scores(path):
    """Read a two-column table (variant ID, numeric score); returns a Series indexed by variant.

    Duplicate variants keep the maximum score.
    """
    df = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], dtype=str, comment="#")
    df.columns = ["variant", "score"]
    df["variant"] = normalize_variant_ids(df["variant"])["variant"].values
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    df = df.dropna()
    return df.groupby("variant")["score"].max()
