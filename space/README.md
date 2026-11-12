---
title: finemap-qc
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: gradio
sdk_version: 6.31.0
app_file: app.py
python_version: "3.12"
short_description: Enrichment-based QC for fine-mapped credible sets
startup_duration_timeout: 30m
---

# finemap-qc

Upload fine-mapped credible sets (e.g. eQTL Catalogue `*.credible_sets.tsv.gz`) and get a report on how
strongly high-PIP variants are enriched for features of causal variants compared with low-PIP variants.

Source: https://github.com/kerimoff/finamapped

Reference annotations (GENCODE, ENCODE cCREs, eQTL Catalogue replicating variants) are downloaded at
start-up from the dataset named in the `FINEMAP_QC_REFS_REPO` variable.
