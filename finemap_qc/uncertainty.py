"""Cluster bootstrap for unique-variant enrichment and median differences.

Connected credible sets are a dependence proxy, not inferred LD blocks. A supplied
locus_id joins sets before shared variants are used to merge remaining clusters.
"""
import numpy as np
import pandas as pd


def variant_clusters(cs):
    """Assign each variant to one connected component, deterministically."""
    parent = {c: c for c in sorted(cs.cs_id.unique())}

    def root(c):
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    for column in (["locus_id"] if "locus_id" in cs else []) + ["variant"]:
        first = {}
        for c, key in zip(cs.cs_id, cs[column]):
            if pd.isna(key):
                raise ValueError(f"missing {column} for clustering")
            if key in first:
                a, b = root(c), root(first[key])
                parent[max(a, b)] = min(a, b)
            else:
                first[key] = c
    pairs = cs.drop_duplicates("variant").set_index("variant")
    return pairs.cs_id.map(root).rename("cluster")


def _weighted_median(values, counts):
    """Median of sorted values with integer multiplicities (no expanded arrays)."""
    cumulative = np.cumsum(counts)
    n = cumulative[-1] if len(cumulative) else 0
    if not n:
        return np.nan
    indices = np.searchsorted(cumulative, [(n - 1) // 2, n // 2], side="right")
    return values[indices].mean()


def cluster_intervals(table, annot, groups, clusters, binary, n_boot=1000, seed=0):
    """Add paired cluster percentile intervals, preserving unique-variant estimates.

    Sample K clusters K times with replacement. Repeated draws retain multiplicity;
    every unique variant belongs to exactly one cluster. Missing annotations are
    excluded per feature. Never silently discard undefined replicates to form a CI.
    """
    if n_boot < 100:
        raise ValueError("bootstrap replicates must be >= 100")
    result = table.copy()
    prefix = "fe" if binary else "median_diff"
    for i, row in result.iterrows():
        low = row["comparison"].removeprefix("high_vs_")
        ids = sorted(groups["high"] | groups[low])
        values = annot.loc[ids, row["feature"]].dropna()
        if not np.isfinite(values.to_numpy(dtype=float)).all():
            raise ValueError("annotations must be finite or missing")
        labels = clusters.reindex(values.index)
        if labels.isna().any():
            raise ValueError("missing cluster assignment")
        codes, unique = pd.factorize(labels, sort=True)
        k = len(unique)
        is_high = values.index.isin(groups["high"])
        hi, lo = codes[is_high], codes[~is_high]
        h, l = values.to_numpy(dtype=float)[is_high], values.to_numpy(dtype=float)[~is_high]
        result.loc[i, "ci_method"] = "cluster_percentile"
        result.loc[i, "n_clusters"] = k
        result.loc[i, "n_clusters_high"] = len(set(hi))
        result.loc[i, "n_clusters_low"] = len(set(lo))
        result.loc[i, "bootstrap_replicates"] = n_boot
        result.loc[i, "bootstrap_seed"] = seed
        result.loc[i, f"{prefix}_ci_low"] = np.nan
        result.loc[i, f"{prefix}_ci_high"] = np.nan
        estimates = []
        status = "insufficient_clusters"
        if k >= 2 and len(set(hi)) >= 2 and len(set(lo)) >= 2:
            rng = np.random.default_rng(seed)
            if binary:
                nh = np.bincount(hi, minlength=k)
                nl = np.bincount(lo, minlength=k)
                sh = np.bincount(hi, weights=h, minlength=k)
                sl = np.bincount(lo, weights=l, minlength=k)
            else:
                ih, il = np.argsort(h), np.argsort(l)
                h, hi, l, lo = h[ih], hi[ih], l[il], lo[il]
            for _ in range(n_boot):
                weights = np.bincount(rng.integers(k, size=k), minlength=k)
                if binary:
                    dh, dl = weights @ nh, weights @ nl
                    a, c = weights @ sh, weights @ sl
                    value = (a / dh) / (c / dl) if dh and dl and c else np.nan
                else:
                    value = _weighted_median(h, weights[hi]) - _weighted_median(l, weights[lo])
                estimates.append(value)
            status = "ok" if np.isfinite(estimates).all() else "undefined_replicates"
            if status == "ok":
                result.loc[i, [f"{prefix}_ci_low", f"{prefix}_ci_high"]] = np.quantile(estimates, [.025, .975])
        result.loc[i, "bootstrap_valid"] = int(np.isfinite(estimates).sum())
        result.loc[i, "ci_status"] = status
    return result
