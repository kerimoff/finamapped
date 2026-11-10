"""Plots and the HTML report."""
import base64
import html
import io

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

COLORS = {"high": "#d62728", "low_all": "#1f77b4", "low_anchored": "#2ca02c"}


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def plot_quantitative(annot, groups, feature, path, log_distance=False):
    """Overlaid density histograms per group with dashed lines at the medians."""
    fig, ax = plt.subplots(figsize=(7, 4))
    data = {g: annot.loc[list(vs), feature].dropna().to_numpy() for g, vs in groups.items()}
    transform = (lambda x: np.log10(x + 1)) if log_distance else (lambda x: x)
    pooled = np.concatenate([transform(v) for v in data.values() if len(v)]) if any(len(v) for v in data.values()) else np.array([0.0])
    bins = np.histogram_bin_edges(pooled, bins=40)
    for g, v in data.items():
        if not len(v):
            continue
        med = np.median(v)
        ax.hist(transform(v), bins=bins, density=True, histtype="step", lw=1.5, color=COLORS.get(g),
                label=f"{g} (n={len(v):,}, median={med:,.3g})")
        ax.axvline(transform(med), color=COLORS.get(g), ls="--", lw=1)
    ax.set_xlabel(f"log10({feature} + 1)" if log_distance else feature)
    ax.set_ylabel("density")
    ax.legend(fontsize=8)
    _save(fig, path)


def plot_binary(table, path):
    """Bar chart of the fraction of variants overlapping each binary feature, per group."""
    feats = list(dict.fromkeys(table["feature"]))
    high = table.groupby("feature")["frac_high"].first()
    lows = table.pivot(index="feature", columns="comparison", values="frac_low")
    series = {"high": high}
    for comp in lows.columns:
        series[comp.replace("high_vs_", "")] = lows[comp]
    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(feats) + 2), 4))
    width = 0.8 / len(series)
    x = np.arange(len(feats))
    for i, (g, s) in enumerate(series.items()):
        ax.bar(x + i * width, s.reindex(feats).to_numpy(), width, label=g, color=COLORS.get(g))
    ax.set_xticks(x + width * (len(series) - 1) / 2)
    ax.set_xticklabels(feats, rotation=30, ha="right")
    ax.set_ylabel("fraction of variants")
    ax.legend(fontsize=8)
    _save(fig, path)


def plot_cs_sizes(per_cs, path):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    sizes = per_cs["size"]
    bins = np.unique(np.logspace(0, np.log10(max(sizes.max(), 2)) + 0.01, 30).astype(int))
    axes[0].hist(sizes, bins=bins, color="grey")
    axes[0].set_xscale("log")
    axes[0].axvline(sizes.median(), color="k", ls="--", label=f"median={sizes.median():g}")
    axes[0].set_xlabel("credible set size (variants)")
    axes[0].set_ylabel("credible sets")
    axes[0].legend()
    counts = per_cs.groupby("trait").size().value_counts().sort_index()
    axes[1].bar(counts.index.astype(str), counts.to_numpy(), color="grey")
    axes[1].set_xlabel("credible sets per molecular trait")
    axes[1].set_ylabel("molecular traits")
    _save(fig, path)


def _img(path):
    with open(path, "rb") as fh:
        return f'<img src="data:image/png;base64,{base64.b64encode(fh.read()).decode()}">'


def _table(df):
    return df.to_html(index=False, float_format=lambda x: f"{x:.4g}", border=0, classes="t")


def write_html(path, title, params, summary, binary, quant, images):
    """images: dict section -> list of PNG paths."""
    css = ("body{font-family:sans-serif;margin:2em;max-width:1200px}"
           "table.t{border-collapse:collapse;font-size:13px}"
           "table.t td,table.t th{padding:3px 8px;border-bottom:1px solid #ddd;text-align:right}"
           "img{max-width:100%;margin:8px 0}")
    parts = [f"<html><head><meta charset='utf-8'><title>{html.escape(title)}</title>"
             f"<style>{css}</style></head><body><h1>{html.escape(title)}</h1>",
             "<h2>Parameters</h2>", _table(pd.DataFrame(list(params.items()), columns=["parameter", "value"])),
             "<h2>Descriptive statistics</h2>",
             _table(summary.reset_index().set_axis(["metric", "value"], axis=1))]
    parts += [_img(p) for p in images.get("descriptive", [])]
    if binary is not None and len(binary):
        parts += ["<h2>Binary features</h2>",
                  "<p>frac_* = fraction of unique variants in the group overlapping the feature; "
                  "fold_enrichment = frac_high / frac_low (95% Wald CI). "
                  "P-values assume independent variants and are therefore anti-conservative "
                  "(low-PIP variants are in LD with each other).</p>", _table(binary)]
        parts += [_img(p) for p in images.get("binary", [])]
    if quant is not None and len(quant):
        parts += ["<h2>Quantitative features</h2>", _table(quant)]
        parts += [_img(p) for p in images.get("quantitative", [])]
    parts.append("</body></html>")
    with open(path, "w") as fh:
        fh.write("\n".join(parts))
