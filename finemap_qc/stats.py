"""Variant groups, descriptive statistics and enrichment tests."""
import numpy as np
import pandas as pd
from scipy import stats


def define_groups(cs, high=0.8, low=0.05, anchor=0.5):
    """Return {group name: set of unique variant IDs}.

    high         : PIP > high in at least one credible set
    low_all      : PIP < low in some credible set
    low_anchored : PIP < low in a credible set whose top variant has PIP > anchor
    Variants that are high-PIP anywhere are removed from both low groups.
    """
    if not (0 <= low < high <= 1 and 0 <= anchor <= 1):
        raise ValueError("thresholds must satisfy 0 <= low < high <= 1 and 0 <= anchor <= 1")
    cs_max = cs.groupby("cs_id")["pip"].transform("max")
    high_v = set(cs.loc[cs["pip"] > high, "variant"])
    low_all = set(cs.loc[cs["pip"] < low, "variant"]) - high_v
    low_anch = set(cs.loc[(cs["pip"] < low) & (cs_max > anchor), "variant"]) - high_v
    return {"high": high_v, "low_all": low_all, "low_anchored": low_anch}


def describe(cs, groups, anchor=0.5):
    """Descriptive statistics. Returns (summary Series, per-credible-set DataFrame)."""
    per_cs = cs.groupby("cs_id").agg(trait=("trait", "first"), size=("variant", "size"),
                                     max_pip=("pip", "max"))
    cs_per_trait = per_cs.groupby("trait").size()
    s = per_cs["size"]
    summary = {
        "n_rows": len(cs),
        "n_unique_variants": cs["variant"].nunique(),
        "n_credible_sets": len(per_cs),
        "n_molecular_traits": len(cs_per_trait),
        "cs_per_trait_mean": cs_per_trait.mean(),
        "cs_per_trait_median": cs_per_trait.median(),
        "cs_per_trait_max": cs_per_trait.max(),
        "frac_traits_with_1_cs": (cs_per_trait == 1).mean(),
        "cs_size_median": s.median(),
        "cs_size_mean": s.mean(),
        "cs_size_q25": s.quantile(0.25),
        "cs_size_q75": s.quantile(0.75),
        "cs_size_max": s.max(),
        "frac_cs_size_1": (s == 1).mean(),
        "frac_cs_size_le_10": (s <= 10).mean(),
        f"frac_cs_with_pip_gt_{anchor}": (per_cs["max_pip"] > anchor).mean(),
        "frac_cs_with_pip_gt_0.8": (per_cs["max_pip"] > 0.8).mean(),
        "frac_cs_with_pip_gt_0.95": (per_cs["max_pip"] > 0.95).mean(),
    }
    for name, vs in groups.items():
        summary[f"n_variants_{name}"] = len(vs)
    return pd.Series(summary, dtype=object), per_cs


def binary_enrichment(annot, groups, features, low_groups=("low_all", "low_anchored")):
    """Fraction of high- vs low-PIP variants carrying each binary feature.

    fold_enrichment = frac_high / frac_low, with a Wald 95% CI on the log scale.
    odds_ratio and p-value from Fisher's exact test.
    """
    rows = []
    high = list(groups["high"])
    for low_name in low_groups:
        low = list(groups[low_name])
        for f in features:
            h = annot.loc[high, f].dropna().astype(bool)
            l = annot.loc[low, f].dropna().astype(bool)
            a, n1, c, n2 = int(h.sum()), len(h), int(l.sum()), len(l)
            row = {"feature": f, "comparison": f"high_vs_{low_name}",
                   "n_high": n1, "k_high": a, "frac_high": a / n1 if n1 else np.nan,
                   "n_low": n2, "k_low": c, "frac_low": c / n2 if n2 else np.nan,
                   "fold_enrichment": np.nan, "fe_ci_low": np.nan, "fe_ci_high": np.nan,
                   "odds_ratio": np.nan, "fisher_p": np.nan}
            if n1 and n2:
                if c:
                    row["fold_enrichment"] = row["frac_high"] / row["frac_low"]
                if a and c:
                    se = np.sqrt(1 / a - 1 / n1 + 1 / c - 1 / n2)
                    row["fe_ci_low"] = row["fold_enrichment"] * np.exp(-1.96 * se)
                    row["fe_ci_high"] = row["fold_enrichment"] * np.exp(1.96 * se)
                row["odds_ratio"], row["fisher_p"] = stats.fisher_exact([[a, n1 - a], [c, n2 - c]])
            rows.append(row)
    return pd.DataFrame(rows)


def quantitative_comparison(annot, groups, features, low_groups=("low_all", "low_anchored")):
    """Medians of quantitative features in high- vs low-PIP variants (+ Mann-Whitney U test).

    Variants with a missing value (e.g. not in the score table) are dropped per feature.
    """
    rows = []
    for low_name in low_groups:
        for f in features:
            h = annot.loc[list(groups["high"]), f].dropna()
            l = annot.loc[list(groups[low_name]), f].dropna()
            row = {"feature": f, "comparison": f"high_vs_{low_name}",
                   "n_high": len(h), "n_low": len(l),
                   "median_high": h.median(), "median_low": l.median(),
                   "mean_high": h.mean(), "mean_low": l.mean()}
            row["median_diff"] = row["median_high"] - row["median_low"]
            row["median_ratio"] = row["median_high"] / row["median_low"] if row["median_low"] else np.nan
            row["mannwhitney_p"] = stats.mannwhitneyu(h, l).pvalue if len(h) and len(l) else np.nan
            rows.append(row)
    return pd.DataFrame(rows)
