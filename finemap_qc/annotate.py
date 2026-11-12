"""Variant annotation: distances to gene features and interval overlaps."""
from functools import lru_cache

import numpy as np
import pandas as pd


def _attr(attrs, key):
    return attrs.str.extract(rf'{key} "([^"]+)"', expand=False)


@lru_cache(maxsize=4)  # the web app reuses parsed reference files between runs
def load_gene_features(gtf_path, canonical_only=False):
    """Extract TSS and splice-site positions of protein-coding transcripts from a GENCODE GTF.

    Returns {"tss": {chrom: sorted positions}, "splice": {chrom: sorted positions}}.
    Splice sites are the first/last exonic base at every exon-intron boundary
    (i.e. exon starts/ends excluding the transcript's own start and end).
    Coordinates are 1-based, as in the GTF.
    """
    gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None, usecols=[0, 2, 3, 4, 6, 8],
                      names=["chrom", "feature", "start", "end", "strand", "attrs"],
                      dtype={"chrom": str})
    gtf = gtf[gtf["feature"].isin(["transcript", "exon"])]
    gtf["chrom"] = gtf["chrom"].str.replace(r"^chr", "", regex=True)
    gtf["tx"] = _attr(gtf["attrs"], "transcript_id")
    tx = gtf[gtf["feature"] == "transcript"]
    keep = (tx["attrs"].str.contains('gene_type "protein_coding"', regex=False)
            & tx["attrs"].str.contains('transcript_type "protein_coding"', regex=False))
    if canonical_only:
        keep &= tx["attrs"].str.contains('tag "Ensembl_canonical"', regex=False)
    tx = tx[keep]
    # Transcript selection controls exon selection; tags need not be repeated on exons.
    gtf = gtf[gtf["tx"].isin(tx["tx"])]
    tss = pd.DataFrame({"chrom": tx["chrom"], "pos": np.where(tx["strand"] == "+", tx["start"], tx["end"])})

    ex = gtf[gtf["feature"] == "exon"]
    ex = ex.drop_duplicates(["tx", "chrom", "start", "end"]).sort_values(["tx", "start", "end"])
    next_start = ex.groupby("tx")["start"].shift(-1)
    # Only adjacent exons separated by an intron define splice boundaries.
    junction = next_start > ex["end"] + 1
    splice = pd.concat([
        pd.DataFrame({"chrom": ex.loc[junction, "chrom"], "pos": ex.loc[junction, "end"]}),
        pd.DataFrame({"chrom": ex.loc[junction, "chrom"], "pos": next_start[junction].astype("int64")}),
    ])

    def by_chrom(df):
        return {c: np.unique(g["pos"].to_numpy()) for c, g in df.groupby("chrom")}

    return {"tss": by_chrom(tss), "splice": by_chrom(splice)}


def nearest_distance(chroms, pos, sites):
    """Absolute distance from each variant to the nearest site on the same chromosome (NaN if none)."""
    chroms = np.asarray(chroms)
    pos = np.asarray(pos, dtype="int64")
    out = np.full(len(pos), np.nan)
    for c in np.unique(chroms):
        s = sites.get(c)
        if s is None or len(s) == 0:
            continue
        m = chroms == c
        p = pos[m]
        i = np.searchsorted(s, p)
        left = s[np.clip(i - 1, 0, len(s) - 1)]
        right = s[np.clip(i, 0, len(s) - 1)]
        out[m] = np.minimum(np.abs(p - left), np.abs(p - right))
    return out


def overlaps_intervals(chroms, pos, bed):
    """True if the variant position (1-based) falls inside any BED interval (0-based, half-open)."""
    chroms = np.asarray(chroms)
    p0 = np.asarray(pos, dtype="int64") - 1
    out = np.zeros(len(p0), dtype=bool)
    for c in np.unique(chroms):
        if c not in bed:
            continue
        starts, max_ends = bed[c]
        m = chroms == c
        i = np.searchsorted(starts, p0[m], side="right") - 1
        hit = np.zeros(m.sum(), dtype=bool)
        ok = i >= 0
        hit[ok] = max_ends[i[ok]] > p0[m][ok]
        out[m] = hit
    return out
