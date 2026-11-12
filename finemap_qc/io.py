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
    # object dtype keeps .str usable when input is empty or no ID splits into 4 parts
    parts = parts.reindex(columns=range(4)).astype(object)
    chrom = parts[0].str.replace(r"^chr", "", regex=True)
    pos = pd.to_numeric(parts[1], errors="coerce")
    ok = (parts.notna().all(axis=1) & np.isfinite(pos) & (pos > 0)
          & (pos % 1 == 0) & (pos < 2**63)
          & parts[1].str.fullmatch(r"[0-9]+", na=False)
          & chrom.ne("") & parts[2].str.fullmatch(r"[ACGTNacgtn]+", na=False)
          & parts[3].str.fullmatch(r"[ACGTNacgtn]+", na=False))
    parts[1] = parts[1].str.lstrip("0")
    parts[2] = parts[2].str.upper()
    parts[3] = parts[3].str.upper()
    variant = ("chr" + chrom + "_" + parts[1] + "_" + parts[2] + "_" + parts[3]).where(ok)
    return pd.DataFrame({"variant": variant, "chrom": chrom.where(ok), "pos": pos.where(ok)})


def read_credible_sets(path, duplicate_policy="error"):
    """Read credible sets.

    Uses columns named cs_id, variant, pip if present (so eQTL Catalogue
    *.credible_sets.tsv.gz files work directly); otherwise the first three
    headerless columns. Named headers must include cs_id, variant and pip.  If a molecular_trait_id column
    exists it is used to group credible sets by trait; otherwise the trait is
    the cs_id with a trailing ``_L<number>`` removed (eQTL Catalogue convention).
    """
    try:
        raw = pd.read_csv(path, sep="\t", header=None, dtype=str)
    except pd.errors.EmptyDataError as exc:
        raise ValueError("credible-set input is empty") from exc
    if raw.shape[1] < 3:
        raise ValueError("credible-set input requires at least three columns")
    first = raw.iloc[0].tolist()
    required = {"cs_id", "variant", "pip"}
    # Recognize named headers, never infer a header from a malformed data PIP.
    has_header = bool(required.intersection(first))
    if has_header:
        if len(set(first)) != len(first):
            raise ValueError("duplicate column names in credible-set header")
        if not required <= set(first):
            raise ValueError("header must contain cs_id, variant and pip")
        raw = raw.iloc[1:].reset_index(drop=True)
        raw.columns = first
        keep = ["cs_id", "variant", "pip"] + [
            c for c in ("molecular_trait_id", "locus_id") if c in raw.columns]
        df = raw[keep].copy()
    else:
        df = raw.iloc[:, :3].copy()
        df.columns = ["cs_id", "variant", "pip"]
    if df.empty:
        raise ValueError("credible-set input has no data rows")
    for col in ("cs_id", "molecular_trait_id", "locus_id"):
        if col in df and (df[col].isna() | df[col].str.strip().eq("")).any():
            raise ValueError(f"{col} must be non-missing and non-empty")
    df["pip"] = pd.to_numeric(df["pip"], errors="raise")
    if not (np.isfinite(df["pip"]) & df["pip"].between(0, 1)).all():
        raise ValueError("PIP must be finite and in [0, 1]")
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
    if duplicate_policy not in ("error", "collapse-identical"):
        raise ValueError("unknown duplicate policy")
    duplicates = df.duplicated(["cs_id", "variant"])
    removed = 0
    if duplicates.any():
        collapsed = df.drop_duplicates()
        conflicting = collapsed.duplicated(["cs_id", "variant"], keep=False)
        if conflicting.any():
            ex = collapsed.loc[conflicting].iloc[0]
            raise ValueError(
                f"duplicate credible-set/variant rows with conflicting values: {conflicting.sum()} rows, "
                f"e.g. {ex['variant']} in {ex['cs_id']} (PIP, trait or locus differ). Each variant must "
                "appear once per credible set.")
        if duplicate_policy == "error":
            ex = df.loc[duplicates].iloc[0]
            raise ValueError(
                f"duplicate credible-set/variant rows: {int(duplicates.sum())} repeated rows, e.g. "
                f"{ex['variant']} in {ex['cs_id']}. eQTL Catalogue files list a variant once per rsID; "
                "the rows are otherwise identical and can be merged with --duplicate-policy collapse-identical "
                "(web app: 'Collapse repeated variant rows').")
        removed = len(df) - len(collapsed)
        df = collapsed.reset_index(drop=True)
    df.attrs["duplicate_rows_removed"] = removed
    for col in ("trait", "locus_id"):
        if col in df and (df.groupby("cs_id")[col].nunique() > 1).any():
            raise ValueError(f"inconsistent {col} IDs within a credible set")
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
