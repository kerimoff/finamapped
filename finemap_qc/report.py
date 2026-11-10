"""Plots and the HTML report.

Every figure is saved as a PNG (for the CLI output folder and the web UI) and as an
SVG next to it.  The report inlines the SVG and swaps the fixed "sentinel" colours
below for CSS variables, so the charts follow the page's light/dark theme.
"""
import datetime as _dt
import functools
import html
import math
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import ticker  # noqa: E402

# Sentinel colours used when drawing. Each maps to a CSS variable in the report.
INK_1, INK_2, AXIS, GRID, SURFACE, NEUTRAL = "#1b1b19", "#5d5c57", "#b9b8b1", "#e6e5e0", "#fcfcfb", "#8f8e88"
COLORS = {"high": "#2a78d6", "low_all": "#eb6834", "low_anchored": "#1baf7a"}
_CSS_VARS = {INK_1: "--ink-1", INK_2: "--ink-2", AXIS: "--axis", GRID: "--grid", SURFACE: "--surface",
             NEUTRAL: "--neutral", **{c: f"--g-{g.replace('_', '-')}" for g, c in COLORS.items()}}

_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 11,
    "text.color": INK_1,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": INK_2,
    "axes.labelsize": 11,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.grid.axis": "y",
    "axes.axisbelow": True,
    "axes.facecolor": "none",
    "grid.color": GRID,
    "grid.linewidth": 0.8,
    "xtick.color": AXIS,
    "ytick.color": AXIS,
    "xtick.labelcolor": INK_2,
    "ytick.labelcolor": INK_2,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.frameon": False,
    "legend.fontsize": 10,
    "legend.labelcolor": INK_1,
    "svg.fonttype": "none",
}


def _styled(fn):
    """Run a plotting function with the report's matplotlib style."""
    @functools.wraps(fn)
    def wrapper(*args, **kw):
        with plt.rc_context(_RC):
            return fn(*args, **kw)
    return wrapper


_subplots = plt.subplots


def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    with plt.rc_context({"svg.hashsalt": os.path.basename(path), "svg.fonttype": "none"}):
        fig.savefig(os.path.splitext(path)[0] + ".svg", transparent=True)
    plt.close(fig)


def _legend(ax, ncol=3):
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.01), ncol=ncol, handlelength=1.4,
              columnspacing=1.4, borderaxespad=0)


@_styled
def plot_quantitative(annot, groups, feature, path, log_distance=False):
    """Overlaid density histograms per group with dashed lines at the medians."""
    fig, ax = _subplots(figsize=(5.4, 3.5))
    data = {g: annot.loc[list(vs), feature].dropna().to_numpy() for g, vs in groups.items()}
    transform = (lambda x: np.log10(x + 1)) if log_distance else (lambda x: x)
    pooled = np.concatenate([transform(v) for v in data.values() if len(v)]) if any(len(v) for v in data.values()) else np.array([0.0])
    bins = np.histogram_bin_edges(pooled, bins=40)
    for g, v in data.items():
        if not len(v):
            continue
        med = np.median(v)
        ax.hist(transform(v), bins=bins, density=True, histtype="step", lw=2, color=COLORS.get(g, NEUTRAL),
                label=f"{g}  n={len(v):,} · median {_num(med)}")
        ax.axvline(transform(med), color=COLORS.get(g, NEUTRAL), ls=(0, (4, 3)), lw=1.2)
    ax.set_xlabel(f"log10({feature} + 1)" if log_distance else feature)
    ax.set_ylabel("density")
    _legend(ax, ncol=1)
    _save(fig, path)


@_styled
def plot_binary(table, path):
    """Grouped bars: fraction of variants overlapping each binary feature, per group."""
    feats = list(dict.fromkeys(table["feature"]))
    high = table.groupby("feature")["frac_high"].first()
    lows = table.pivot(index="feature", columns="comparison", values="frac_low")
    series = {"high": high}
    for comp in lows.columns:
        series[comp.replace("high_vs_", "")] = lows[comp]
    fig, ax = _subplots(figsize=(max(5.4, 1.1 * len(feats) + 2.4), 3.5))
    width = 0.8 / len(series)
    x = np.arange(len(feats))
    for i, (g, s) in enumerate(series.items()):
        ax.bar(x + i * width, s.reindex(feats).to_numpy(), width, label=g, color=COLORS.get(g, NEUTRAL),
               edgecolor=SURFACE, linewidth=1.5)
    ax.set_xticks(x + width * (len(series) - 1) / 2)
    ax.set_xticklabels(feats, rotation=25 if len(feats) > 4 else 0, ha="right" if len(feats) > 4 else "center")
    ax.tick_params(axis="x", length=0)
    ax.yaxis.set_major_formatter(ticker.PercentFormatter(1, decimals=0))
    ax.set_ylabel("variants overlapping feature")
    _legend(ax)
    _save(fig, path)


@_styled
def plot_forest(table, path):
    """Fold enrichment (high vs each low group) with 95% CIs, on a log scale."""
    t = table[np.isfinite(table[["fold_enrichment", "fe_ci_low", "fe_ci_high"]]).all(axis=1)
              & (table["fe_ci_low"] > 0)]
    feats = list(dict.fromkeys(table["feature"]))
    comps = list(dict.fromkeys(table["comparison"]))
    fig, ax = _subplots(figsize=(5.4, max(3.5, 0.8 * len(feats) + 1.6)))
    offsets = np.linspace(-0.16, 0.16, len(comps)) if len(comps) > 1 else [0.0]
    for off, comp in zip(offsets, comps):
        g = comp.replace("high_vs_", "")
        r = t[t["comparison"] == comp].set_index("feature").reindex(feats)
        y = np.arange(len(feats))[::-1] + off
        fe = r["fold_enrichment"].to_numpy()
        err = np.vstack([fe - r["fe_ci_low"].to_numpy(), r["fe_ci_high"].to_numpy() - fe])
        ax.errorbar(fe, y, xerr=err, fmt="o", ms=6.5, lw=2, capsize=0, color=COLORS.get(g, NEUTRAL),
                    mec=SURFACE, mew=1.5, label=f"high vs {g}")
    ax.axvline(1, color=INK_2, lw=1, ls=(0, (4, 3)))
    ax.set_xscale("log")
    lo = np.nanmin(t["fe_ci_low"]) if len(t) else 0.5
    hi = np.nanmax(t["fe_ci_high"]) if len(t) else 2
    ax.set_xlim(min(lo, 1) / 1.4, max(hi, 1) * 1.4)
    ax.xaxis.set_major_locator(ticker.LogLocator(subs=(1, 2, 5)))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(ticker.NullLocator())
    ax.set_yticks(np.arange(len(feats))[::-1])
    ax.set_yticklabels(feats)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.6, len(feats) - 0.4)
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("fold enrichment (95% CI, log scale)")
    _legend(ax, ncol=2)
    _save(fig, path)


@_styled
def plot_cs_sizes(per_cs, path):
    fig, axes = _subplots(1, 2, figsize=(10, 3.6), gridspec_kw={"width_ratios": [1.6, 1]})
    sizes = per_cs["size"]
    bins = np.unique(np.logspace(0, np.log10(max(sizes.max(), 2)) + 0.01, 30).astype(int))
    axes[0].hist(sizes, bins=bins, color=NEUTRAL, edgecolor=SURFACE, linewidth=1)
    axes[0].set_xscale("log")
    axes[0].xaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    axes[0].axvline(sizes.median(), color=INK_1, ls=(0, (4, 3)), lw=1.2, label=f"median {sizes.median():g}")
    axes[0].set_xlabel("credible set size (variants)")
    axes[0].set_ylabel("credible sets")
    _legend(axes[0])
    counts = per_cs.groupby("trait").size().value_counts().sort_index()
    axes[1].bar(counts.index.astype(str), counts.to_numpy(), color=NEUTRAL, edgecolor=SURFACE, linewidth=1.5,
                width=0.7)
    axes[1].tick_params(axis="x", length=0)
    axes[1].set_xlabel("credible sets per molecular trait")
    axes[1].set_ylabel("molecular traits")
    _save(fig, path)


# --------------------------------------------------------------------------- HTML

_CAPTIONS = {
    "cs_sizes": ("Credible set sizes", "Left: number of variants per credible set (log scale). "
                 "Right: how many credible sets each molecular trait has."),
    "forest": ("Fold enrichment", "Enrichment of each feature among high-PIP variants relative to each "
               "low-PIP group. Points right of the dashed line at 1 mean the feature is more common "
               "among high-PIP variants."),
    "binary_features": ("Feature overlap by group", "Share of unique variants in each group that "
                        "overlap the feature."),
}


def _caption(path):
    stem = os.path.splitext(os.path.basename(path))[0]
    if stem in _CAPTIONS:
        return _CAPTIONS[stem]
    if stem.startswith("hist_"):
        f = stem[5:]
        return f, f"Distribution of {f} per group; dashed lines mark each group's median."
    return stem, ""


def _svg(path):
    svg_path = os.path.splitext(path)[0] + ".svg"
    if not os.path.exists(svg_path):
        import base64
        with open(path, "rb") as fh:
            return f'<img alt="" src="data:image/png;base64,{base64.b64encode(fh.read()).decode()}">'
    with open(svg_path) as fh:
        s = fh.read()
    s = s[s.index("<svg"):]
    s = re.sub(r'<svg([^>]*?)\swidth="[^"]*"\sheight="[^"]*"', r'<svg\1', s, count=1)
    s = s.replace("<svg", '<svg role="img" preserveAspectRatio="xMidYMid meet"', 1)

    def swap(m):
        var = _CSS_VARS.get(m.group(0).lower())
        return f"var({var})" if var else m.group(0)
    return re.sub(r"#[0-9a-fA-F]{6}\b", swap, s)


def _figure(path):
    title, text = _caption(path)
    return (f'<figure class="card fig"><figcaption><h3>{html.escape(title)}</h3>'
            f'<p>{html.escape(text)}</p></figcaption>{_svg(path)}</figure>')


def _isnan(x):
    try:
        return x is None or (isinstance(x, float) and math.isnan(x)) or pd.isna(x)
    except (TypeError, ValueError):
        return False


def _int(x):
    return "—" if _isnan(x) else f"{int(round(float(x))):,}"


def _num(x, sig=3):
    if _isnan(x):
        return "—"
    x = float(x)
    if math.isinf(x):
        return "∞" if x > 0 else "−∞"
    if abs(x) >= 1000:
        return f"{x:,.0f}".replace("-", "−")
    if x == int(x) and abs(x) < 1000:
        return f"{int(x)}".replace("-", "−")
    return f"{x:.{sig}g}".replace("-", "−")


def _pct(x, d=1):
    return "—" if _isnan(x) else f"{100 * float(x):.{d}f}%"


def _p(p):
    if _isnan(p):
        return "—"
    p = float(p)
    if p == 0:
        return "&lt; 10<sup>−300</sup>"
    if p >= 0.001:
        return f"{p:.3f}"
    exp = math.floor(math.log10(p))
    return f"{p / 10 ** exp:.1f} × 10<sup>−{-exp}</sup>"


def _chip(group):
    cls = html.escape(group.replace("_", "-"))
    return f'<span class="chip"><i class="sw g-{cls}"></i>{html.escape(group)}</span>'


def _comparison(comp):
    return f'<span class="vs">high vs</span> {_chip(comp.replace("high_vs_", ""))}'


def _sig(p):
    return "" if _isnan(p) or float(p) >= 0.05 else ' class="sig"'


def _feature_rows(df):
    """Yield (row_html_cells, feature, is_first_of_feature) so feature names span rows."""
    order = {f: i for i, f in enumerate(dict.fromkeys(df["feature"]))}
    df = df.iloc[np.argsort(df["feature"].map(order).to_numpy(), kind="stable")]
    counts = df["feature"].value_counts()
    seen = set()
    for _, r in df.iterrows():
        first = r["feature"] not in seen
        seen.add(r["feature"])
        cell = (f'<th scope="rowgroup" rowspan="{counts[r["feature"]]}" class="feat">{html.escape(str(r["feature"]))}</th>'
                if first else "")
        yield r, cell, first


def _binary_table(df):
    head = ("<tr><th>Feature</th><th class='l'>Comparison</th><th>High-PIP overlap</th>"
            "<th>Low-PIP overlap</th><th>Fold enrichment <span class='u'>95% CI</span></th>"
            "<th>Odds ratio</th><th>Fisher p</th></tr>")
    rows = []
    for r, feat, first in _feature_rows(df):
        rows.append(
            f"<tr class='{'first' if first else ''}'>{feat}<td class='l'>{_comparison(r['comparison'])}</td>"
            f"<td>{_pct(r['frac_high'])} <span class='u'>{_int(r['k_high'])}/{_int(r['n_high'])}</span></td>"
            f"<td>{_pct(r['frac_low'])} <span class='u'>{_int(r['k_low'])}/{_int(r['n_low'])}</span></td>"
            f"<td><b>{_num(r['fold_enrichment'])}</b> <span class='u'>{_num(r['fe_ci_low'])}–{_num(r['fe_ci_high'])}</span></td>"
            f"<td>{_num(r['odds_ratio'])}</td><td{_sig(r['fisher_p'])}>{_p(r['fisher_p'])}</td></tr>")
    return f"<div class='tw'><table class='t'><thead>{head}</thead><tbody>{''.join(rows)}</tbody></table></div>"


def _quant_table(df):
    head = ("<tr><th>Feature</th><th class='l'>Comparison</th><th>Median high <span class='u'>n</span></th>"
            "<th>Median low <span class='u'>n</span></th><th>Median ratio</th><th>Mean high</th>"
            "<th>Mean low</th><th>Mann–Whitney p</th></tr>")
    rows = []
    for r, feat, first in _feature_rows(df):
        rows.append(
            f"<tr class='{'first' if first else ''}'>{feat}<td class='l'>{_comparison(r['comparison'])}</td>"
            f"<td><b>{_num(r['median_high'])}</b> <span class='u'>{_int(r['n_high'])}</span></td>"
            f"<td><b>{_num(r['median_low'])}</b> <span class='u'>{_int(r['n_low'])}</span></td>"
            f"<td>{_num(r['median_ratio'])}</td><td>{_num(r['mean_high'])}</td><td>{_num(r['mean_low'])}</td>"
            f"<td{_sig(r['mannwhitney_p'])}>{_p(r['mannwhitney_p'])}</td></tr>")
    return f"<div class='tw'><table class='t'><thead>{head}</thead><tbody>{''.join(rows)}</tbody></table></div>"


_SUMMARY_SECTIONS = [
    ("Input", ["n_rows", "n_unique_variants", "n_credible_sets", "n_molecular_traits"]),
    ("Credible set size", ["cs_size_median", "cs_size_mean", "cs_size_q25", "cs_size_q75", "cs_size_max",
                           "frac_cs_size_1", "frac_cs_size_le_10"]),
    ("Resolution", ["frac_cs_with_pip_gt_*"]),
    ("Credible sets per trait", ["cs_per_trait_mean", "cs_per_trait_median", "cs_per_trait_max",
                                 "frac_traits_with_1_cs"]),
]
_LABELS = {
    "n_rows": "Rows in input", "n_unique_variants": "Unique variants", "n_credible_sets": "Credible sets",
    "n_molecular_traits": "Molecular traits", "cs_size_median": "Median", "cs_size_mean": "Mean",
    "cs_size_q25": "25th percentile", "cs_size_q75": "75th percentile", "cs_size_max": "Maximum",
    "frac_cs_size_1": "Single-variant sets", "frac_cs_size_le_10": "Sets with ≤ 10 variants",
    "cs_per_trait_mean": "Mean", "cs_per_trait_median": "Median", "cs_per_trait_max": "Maximum",
    "frac_traits_with_1_cs": "Traits with one set",
}


def _summary_value(key, v):
    if key.startswith("frac_"):
        return _pct(v)
    if key.startswith("n_") or key.endswith("_max"):
        return _int(v)
    return _num(v)


def _summary_label(key):
    if key in _LABELS:
        return html.escape(_LABELS[key])
    m = re.match(r"frac_cs_with_pip_gt_(.+)", key)
    if m:
        return f"Top PIP &gt; {html.escape(m.group(1))}"
    if key.startswith("n_variants_"):
        return _chip(key[len("n_variants_"):])
    return html.escape(key)


def _summary_cards(summary):
    keys = [k for k in summary.index if not k.startswith("n_variants_")]  # shown under "Variant groups"
    used, out = set(), []
    for title, patterns in _SUMMARY_SECTIONS + [("Other", ["*"])]:
        sel = [k for p in patterns for k in keys
               if k not in used and (k == p or (p.endswith("*") and k.startswith(p[:-1])))]
        if not sel:
            continue
        used.update(sel)
        rows = "".join(f"<div class='kv'><dt>{_summary_label(k)}</dt><dd>{_summary_value(k, summary[k])}</dd></div>"
                       for k in sel)
        out.append(f"<section class='card stat'><h3>{html.escape(title)}</h3><dl>{rows}</dl></section>")
    return f"<div class='stats'>{''.join(out)}</div>"


def _kpis(summary):
    g = summary.get
    tiles = [("Credible sets", _int(g("n_credible_sets")), f"{_int(g('n_molecular_traits'))} molecular traits"),
             ("Unique variants", _int(g("n_unique_variants")), f"{_int(g('n_rows'))} rows"),
             ("Median set size", _num(g("cs_size_median")),
              f"IQR {_num(g('cs_size_q25'))}–{_num(g('cs_size_q75'))} variants"),
             ("Sets with top PIP > 0.8", _pct(g("frac_cs_with_pip_gt_0.8")),
              f"{_pct(g('frac_cs_with_pip_gt_0.95'))} above 0.95"),
             ("Single-variant sets", _pct(g("frac_cs_size_1")), f"{_pct(g('frac_cs_size_le_10'))} have ≤ 10 variants")]
    return "<div class='kpis'>" + "".join(
        f"<div class='card kpi'><div class='kl'>{html.escape(a)}</div><div class='kn'>{b}</div>"
        f"<div class='ks'>{html.escape(c)}</div></div>" for a, b, c in tiles) + "</div>"


def _groups(params, summary):
    defs = {"high": ("PIP", params.get("high_pip", "")),
            "low_all": ("PIP", params.get("low_pip", "")),
            "low_anchored": ("PIP", f"{params.get('low_pip', '')}, in a set whose top PIP is "
                             f"{params.get('anchor_pip (low_anchored)', '')}")}
    items = []
    for g, (lab, rule) in defs.items():
        n = summary.get(f"n_variants_{g}")
        items.append(f"<li>{_chip(g)}<span class='gd'>{lab} {html.escape(str(rule))}</span>"
                     f"<span class='gn'>{_int(n)} variants</span></li>")
    return f"<ul class='groups'>{''.join(items)}</ul>"


_CSS = """
:root{color-scheme:light;
--bg:#f4f3ef;--surface:#fcfcfb;--ink-1:#1b1b19;--ink-2:#5d5c57;--ink-3:#8a8983;--line:#e6e5e0;
--axis:#b9b8b1;--grid:#e6e5e0;--neutral:#8f8e88;--accent:#2a78d6;--sig-bg:#e6f0fc;--sig:#1c5cab;
--g-high:#2a78d6;--g-low-all:#eb6834;--g-low-anchored:#1baf7a;--shadow:0 1px 2px rgba(20,20,15,.05),0 2px 10px rgba(20,20,15,.04)}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){color-scheme:dark;
--bg:#121211;--surface:#1c1c1a;--ink-1:#f0efec;--ink-2:#c3c2b7;--ink-3:#8f8e86;--line:#2e2e2b;
--axis:#55544f;--grid:#2e2e2b;--neutral:#8a8983;--accent:#3987e5;--sig-bg:#17304f;--sig:#9ec5f4;
--g-high:#3987e5;--g-low-all:#d95926;--g-low-anchored:#199e70;--shadow:none}}
:root[data-theme="dark"]{color-scheme:dark;
--bg:#121211;--surface:#1c1c1a;--ink-1:#f0efec;--ink-2:#c3c2b7;--ink-3:#8f8e86;--line:#2e2e2b;
--axis:#55544f;--grid:#2e2e2b;--neutral:#8a8983;--accent:#3987e5;--sig-bg:#17304f;--sig:#9ec5f4;
--g-high:#3987e5;--g-low-all:#d95926;--g-low-anchored:#199e70;--shadow:none}
*{box-sizing:border-box}
html{scroll-behavior:smooth;scroll-padding-top:72px}
body{margin:0;background:var(--bg);color:var(--ink-1);
font:15px/1.55 ui-sans-serif,-apple-system,BlinkMacSystemFont,"Segoe UI",Inter,"Helvetica Neue",Arial,sans-serif;
-webkit-font-smoothing:antialiased}
.top{position:sticky;top:0;z-index:5;background:color-mix(in srgb,var(--bg) 86%,transparent);
backdrop-filter:saturate(1.4) blur(10px);-webkit-backdrop-filter:saturate(1.4) blur(10px);border-bottom:1px solid var(--line)}
.top .in{display:flex;align-items:center;gap:20px;height:56px}
.brand{font-weight:650;letter-spacing:-.01em;white-space:nowrap;display:flex;align-items:center;gap:9px}
.brand .dot{width:10px;height:10px;border-radius:3px;background:var(--accent);box-shadow:14px 0 0 var(--g-low-all),28px 0 0 var(--g-low-anchored);margin-right:30px}
nav{display:flex;gap:4px;overflow-x:auto;scrollbar-width:none;margin-left:auto}
nav a{color:var(--ink-2);text-decoration:none;font-size:13.5px;padding:6px 10px;border-radius:8px;white-space:nowrap}
nav a:hover{background:var(--surface);color:var(--ink-1)}
.in{max-width:1180px;margin:0 auto;padding:0 24px}
header.hero{padding:44px 0 8px}
.eyebrow{font-size:12.5px;font-weight:600;letter-spacing:.08em;text-transform:uppercase;color:var(--ink-3)}
h1{font-size:34px;line-height:1.15;letter-spacing:-.025em;margin:6px 0 10px;font-weight:700}
.meta{color:var(--ink-2);font-size:14px;display:flex;flex-wrap:wrap;gap:6px 18px}
.meta code{font:12.5px ui-monospace,SFMono-Regular,Menlo,monospace;background:var(--surface);border:1px solid var(--line);padding:1px 6px;border-radius:6px;word-break:break-all}
section.block{padding:36px 0 4px}
h2{font-size:22px;letter-spacing:-.015em;margin:0 0 6px;font-weight:650}
.lede{color:var(--ink-2);margin:0 0 18px;max-width:760px}
h3{font-size:14.5px;margin:0;font-weight:620}
.card{background:var(--surface);border:1px solid var(--line);border-radius:14px;box-shadow:var(--shadow)}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:12px;margin:26px 0 0}
.kpi{padding:16px 18px}
.kl{font-size:12.5px;color:var(--ink-2);font-weight:550}
.kn{font-size:30px;font-weight:680;letter-spacing:-.02em;margin:4px 0 2px;font-variant-numeric:tabular-nums}
.ks{font-size:12.5px;color:var(--ink-3)}
.groups{list-style:none;padding:0;margin:14px 0 0;display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:12px}
.groups li{background:var(--surface);border:1px solid var(--line);border-radius:12px;padding:12px 16px;display:flex;flex-direction:column;gap:3px}
.gd{font-size:13.5px;color:var(--ink-2)}.gn{font-size:12.5px;color:var(--ink-3);font-variant-numeric:tabular-nums}
.chip{display:inline-flex;align-items:center;gap:7px;font:600 13px ui-monospace,SFMono-Regular,Menlo,monospace;white-space:nowrap}
.sw{display:inline-block;width:10px;height:10px;border-radius:3px;flex:none}
.g-high{background:var(--g-high)}.g-low-all{background:var(--g-low-all)}.g-low-anchored{background:var(--g-low-anchored)}
.vs{color:var(--ink-3);font-size:12.5px;margin-right:2px}
.stats{display:grid;grid-template-columns:repeat(auto-fit,minmax(250px,1fr));gap:12px}
.stat{padding:16px 18px}
.stat dl{margin:10px 0 0}
.kv{display:flex;justify-content:space-between;gap:12px;padding:6px 0;border-top:1px solid var(--line);font-size:14px}
.kv dt{color:var(--ink-2)}.kv dd{margin:0;font-variant-numeric:tabular-nums;font-weight:560}
.figs{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,520px),1fr));gap:14px;margin-top:14px}
.fig{margin:0;padding:18px 18px 10px}
.fig.wide{grid-column:1/-1}
.fig p{margin:3px 0 8px;color:var(--ink-2);font-size:13px}
.fig svg,.fig img{display:block;width:100%;height:auto}
.fig svg text{font-family:inherit!important}
.tw{overflow-x:auto;background:var(--surface);border:1px solid var(--line);border-radius:14px;box-shadow:var(--shadow)}
table.t{border-collapse:collapse;width:100%;font-size:13.5px;font-variant-numeric:tabular-nums}
.t th,.t td{padding:9px 14px;text-align:right;white-space:nowrap;vertical-align:middle}
.t thead th{font-size:12px;font-weight:600;color:var(--ink-2);text-transform:none;border-bottom:1px solid var(--line);background:color-mix(in srgb,var(--bg) 45%,var(--surface))}
.t tbody tr.first{border-top:1px solid var(--line)}
.t tbody tr:first-child{border-top:none}
.t tbody tr:hover td{background:color-mix(in srgb,var(--accent) 5%,transparent)}
.t .l{text-align:left}
.t th.feat{text-align:left;font-weight:620;vertical-align:top;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:13px}
.u{color:var(--ink-3);font-size:12px;font-weight:400;margin-left:3px}
td.sig{color:var(--sig);font-weight:620}
td.sig::after{content:"";display:inline-block;width:6px;height:6px;border-radius:50%;background:var(--sig);margin-left:7px;vertical-align:2px}
.note{display:flex;gap:12px;padding:12px 16px;margin:14px 0 0;border-radius:12px;background:var(--surface);border:1px solid var(--line);border-left:3px solid var(--g-low-all);font-size:13.5px;color:var(--ink-2)}
.legend-sig{font-size:12.5px;color:var(--ink-3);margin:8px 2px 0}
.legend-sig i{display:inline-block;width:6px;height:6px;border-radius:50%;background:var(--sig);margin-right:6px;vertical-align:1px}
.params{padding:6px 18px}
.params .kv dt{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:13px}
.params .kv dd{font:13px ui-monospace,SFMono-Regular,Menlo,monospace;font-weight:400;word-break:break-all;text-align:right}
.params .kv:first-child{border-top:none}
footer{color:var(--ink-3);font-size:12.5px;padding:40px 0 48px}
@media (max-width:640px){.in{padding:0 16px}h1{font-size:27px}.kn{font-size:25px}nav{display:none}header.hero{padding-top:28px}
.params .kv{flex-direction:column;gap:0}.params .kv dd{text-align:left}}
@media print{.top{position:static}.card,.tw{box-shadow:none}}
"""


def write_html(path, title, params, summary, binary, quant, images):
    """images: dict section -> list of PNG paths (an SVG with the same stem is inlined if present)."""
    nav = [("overview", "Overview"), ("credible-sets", "Credible sets")]
    has_bin = binary is not None and len(binary)
    has_q = quant is not None and len(quant)
    if has_bin:
        nav.append(("binary", "Binary features"))
    if has_q:
        nav.append(("quantitative", "Quantitative features"))
    nav.append(("parameters", "Parameters"))
    when = _dt.datetime.now().strftime("%d %b %Y, %H:%M")
    src = params.get("credible_sets", "")

    parts = [
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>{html.escape(title)}</title><style>{_CSS}</style></head><body>",
        "<div class='top'><div class='in'><div class='brand'><span class='dot'></span>finemap-qc</div><nav>"
        + "".join(f"<a href='#{a}'>{html.escape(b)}</a>" for a, b in nav) + "</nav></div></div>",
        "<main class='in'>",
        f"<header class='hero' id='overview'><div class='eyebrow'>Quality-control report</div>"
        f"<h1>{html.escape(title)}</h1><div class='meta'><span>Credible sets <code>{html.escape(str(src))}</code></span>"
        f"<span>Generated {html.escape(when)}</span></div>{_kpis(summary)}"
        f"<h3 style='margin-top:26px'>Variant groups</h3>{_groups(params, summary)}</header>",
        "<section class='block' id='credible-sets'><h2>Credible sets</h2>"
        "<p class='lede'>Size and resolution of the fine-mapped credible sets. Smaller sets and higher top "
        "PIPs mean better-resolved signals.</p>",
        _summary_cards(summary),
        "<div class='figs'>" + "".join(_figure(p).replace("card fig", "card fig wide", 1)
                                       for p in images.get("descriptive", [])) + "</div></section>",
    ]
    if has_bin:
        parts += ["<section class='block' id='binary'><h2>Binary features</h2>"
                  "<p class='lede'>Does a feature (open chromatin, missense, replication…) occur more often among "
                  "high-PIP variants than among low-PIP variants? Well-calibrated fine-mapping should show "
                  "clear enrichment for functional features.</p>",
                  _binary_table(binary),
                  "<p class='legend-sig'><i></i>p &lt; 0.05. Overlap shows the percentage with the count of "
                  "overlapping / total unique variants. Fold enrichment = high-PIP overlap ÷ low-PIP overlap "
                  "(95% Wald CI).</p>",
                  "<div class='note'>P-values treat variants as independent, so they are anti-conservative: "
                  "low-PIP variants in the same credible set are in LD with each other. Read them as a "
                  "ranking, not a calibrated test.</div>",
                  "<div class='figs'>" + "".join(_figure(p) for p in images.get("binary", [])) + "</div></section>"]
    if has_q:
        parts += ["<section class='block' id='quantitative'><h2>Quantitative features</h2>"
                  "<p class='lede'>How continuous annotations differ between high- and low-PIP variants. "
                  "Causal variants are expected to sit closer to genes and splice sites and to score higher "
                  "on functional predictors.</p>",
                  _quant_table(quant),
                  "<p class='legend-sig'><i></i>p &lt; 0.05 (two-sided Mann–Whitney U). Median ratio = "
                  "median high ÷ median low.</p>",
                  "<div class='figs'>" + "".join(_figure(p) for p in images.get("quantitative", [])) + "</div></section>"]
    prow = "".join(f"<div class='kv'><dt>{html.escape(str(k))}</dt><dd>{html.escape(str(v))}</dd></div>"
                   for k, v in params.items())
    parts += [f"<section class='block' id='parameters'><h2>Parameters</h2><div class='card params'><dl>{prow}</dl></div></section>",
              f"<footer>Generated by finemap-qc on {html.escape(when)}.</footer></main></body></html>"]
    with open(path, "w") as fh:
        fh.write("\n".join(parts))
