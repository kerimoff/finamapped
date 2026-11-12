"""Web UI for finemap-qc.  Run with:  streamlit run app.py

Reference annotations (GTF, BED, variant sets, scores) live on the server and are
listed in a JSON file: $FINEMAP_QC_REFS, default examples/references.json.
Users upload only their credible sets (plus, optionally, extra annotation files).
"""
import io
import json
import os
import tempfile
import zipfile
from pathlib import Path

import pandas as pd
import streamlit as st

from finemap_qc import cli

ROOT = Path(__file__).resolve().parent
REFS = Path(os.environ.get("FINEMAP_QC_REFS", ROOT / "examples" / "references.json"))
EXAMPLE = ROOT / "examples" / "credible_sets.tsv.gz"


def reference_args(path):
    """Turn a references.json file into CLI arguments (paths are relative to the JSON file)."""
    cfg = json.loads(path.read_text())
    resolve = lambda p: str((path.parent / p).resolve())  # noqa: E731
    args = ["--gtf", resolve(cfg["gtf"])] if cfg.get("gtf") else []
    for kind in ("bed", "variant_set", "score"):
        for name, p in cfg.get(kind, {}).items():
            args += [f"--{kind.replace('_', '-')}", f"{name}={resolve(p)}"]
    return args, cfg


def save_upload(upload, folder, stem):
    """Write an uploaded file to disk, keeping its extension (so .gz is detected)."""
    path = Path(folder) / (stem + "".join(Path(upload.name).suffixes))
    path.write_bytes(upload.getvalue())
    return path


def zip_folder(folder):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(Path(folder).iterdir()):
            z.write(f, f.name)
    return buf.getvalue()


st.set_page_config(page_title="Fine-mapping QC", layout="wide")
st.title("Fine-mapping QC")
st.write("Upload credible sets (TSV with `cs_id`, `variant` as `chr_pos_ref_alt` on GRCh38, and `pip`; "
         "eQTL Catalogue `*.credible_sets.tsv.gz` files work as they are). The report shows how strongly "
         "high-PIP variants are enriched for features of causal variants, compared with low-PIP variants.")

# ---- sidebar: parameters and annotations -------------------------------------------------------
with st.sidebar:
    st.header("Thresholds")
    high = st.number_input("High PIP >", 0.0, 1.0, 0.8, 0.05)
    low = st.number_input("Low PIP <", 0.0, 1.0, 0.05, 0.01)
    anchor = st.number_input("low_anchored: credible set top PIP >", 0.0, 1.0, 0.5, 0.05)

    st.header("Annotations")
    ref_cli_args, cfg = ([], {})
    if REFS.exists():
        ref_cli_args, cfg = reference_args(REFS)
        use_refs = st.checkbox("Use server reference annotations", value=True)
        st.caption(f"From `{REFS}`")
        st.json(cfg, expanded=False)
        if not use_refs:
            ref_cli_args = []
    else:
        st.warning(f"No reference file found at {REFS}")
    collapse_duplicates = st.checkbox("Collapse identical credible-set/variant rows (Catalogue rsID aliases)", value=True)
    canonical = st.checkbox("GTF: canonical transcripts only", value=False)

    with st.expander("Add your own annotation files"):
        st.caption("The file name (without extension) becomes the feature name.")
        extra_bed = st.file_uploader("BED interval files", accept_multiple_files=True)
        extra_sets = st.file_uploader("Variant ID lists", accept_multiple_files=True)
        extra_scores = st.file_uploader("Score tables (variant, score)", accept_multiple_files=True)

# ---- input ---------------------------------------------------------------------------------------
upload = st.file_uploader("Credible sets file", type=["tsv", "gz", "txt"])
c1, c2 = st.columns([1, 6])
run_upload = c1.button("Run QC", type="primary", disabled=upload is None)
run_example = c2.button(f"Run on example data ({EXAMPLE.name}, simulated)")

if run_upload or run_example:
    work = Path(tempfile.mkdtemp(prefix="finemap_qc_"))
    cs_path = save_upload(upload, work, "credible_sets") if run_upload else EXAMPLE
    argv = ["--credible-sets", str(cs_path), "--out", str(work / "out"),
            "--high-pip", str(high), "--low-pip", str(low), "--anchor-pip", str(anchor)] + ref_cli_args
    if collapse_duplicates:
        argv += ["--duplicate-policy", "collapse-identical"]
    if canonical:
        argv.append("--canonical-only")
    for flag, files in (("--bed", extra_bed), ("--variant-set", extra_sets), ("--score", extra_scores)):
        for f in files or []:
            name = Path(f.name).name.split(".")[0]
            argv += [flag, f"{name}={save_upload(f, work, name + flag)}"]
    with st.spinner("Annotating variants and computing enrichments..."):
        try:
            cli.main(argv)
        except (Exception, SystemExit) as e:  # show input errors in the page instead of crashing
            st.error(f"Run failed: {e}")
            st.stop()
    st.session_state["out"] = str(work / "out")
    st.session_state["input_name"] = upload.name if run_upload else EXAMPLE.name

# ---- results -------------------------------------------------------------------------------------
if "out" in st.session_state:
    out = Path(st.session_state["out"])
    st.subheader(f"Results: {st.session_state['input_name']}")
    summary = pd.read_csv(out / "descriptive_summary.tsv", sep="\t", index_col=0)["value"]
    binary = pd.read_csv(out / "binary_enrichment.tsv", sep="\t")
    quant = pd.read_csv(out / "quantitative_comparison.tsv", sep="\t")

    tab_sum, tab_bin, tab_q, tab_dl = st.tabs(["Credible sets", "Binary features", "Quantitative features",
                                                "Download"])
    with tab_sum:
        cols = st.columns(5)
        cols[0].metric("Credible sets", f"{int(summary['n_credible_sets']):,}")
        cols[1].metric("Molecular traits", f"{int(summary['n_molecular_traits']):,}")
        cols[2].metric("Median CS size", f"{float(summary['cs_size_median']):g}")
        pip_key = next(k for k in summary.index if k.startswith("frac_cs_with_pip_gt_"))
        cols[3].metric(f"CS with PIP > {pip_key.rsplit('_', 1)[1]}", f"{float(summary[pip_key]):.1%}")
        cols[4].metric("High-PIP variants", f"{int(summary['n_variants_high']):,}")
        st.image(str(out / "cs_sizes.png"))
        st.dataframe(summary, width="content")
    with tab_bin:
        if binary.empty:
            st.info("No binary features. Add BED files or variant lists.")
        else:
            st.caption("fold_enrichment = frac_high / frac_low. P-values assume independent variants, "
                       "which is wrong under LD, so they are too small; compare fold enrichments.")
            st.dataframe(binary, hide_index=True)
            st.image(str(out / "binary_features.png"))
    with tab_q:
        if quant.empty:
            st.info("No quantitative features. Add a GTF or score tables.")
        else:
            st.dataframe(quant, hide_index=True)
            for f in quant["feature"].unique():
                st.image(str(out / f"hist_{f}.png"))
    with tab_dl:
        st.download_button("HTML report", (out / "report.html").read_bytes(), "report.html", "text/html")
        st.download_button("All outputs (zip)", zip_folder(out), "finemap_qc_results.zip", "application/zip")
