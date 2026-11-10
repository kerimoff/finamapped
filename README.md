# finemap-qc

Quality metrics for fine-mapped credible sets in eQTL Catalogue format. The tool asks whether
high-PIP variants (PIP > 0.8) are enriched for features that true causal variants tend to have,
compared with low-PIP variants.

## Quick start

```bash
git clone https://github.com/kerimoff/finamapped.git && cd finamapped
git checkout claude/gallant-knuth-tzv77n
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[app,test]"
streamlit run app.py                                     # opens http://localhost:8501 (local only)
```

In the web app, upload a credible sets file, or click **Run on example data**. The results appear in tabs, and you can download the HTML report or a zip of all outputs.

The same analysis is also available from the command line:

```bash
finemap-qc --credible-sets examples/credible_sets.tsv.gz --gtf examples/genes.gtf.gz \
  --bed open_chromatin=examples/open_chromatin.bed --variant-set missense=examples/missense.txt \
  --score spliceai=examples/spliceai.tsv --out qc_out
```

### Example data

`examples/` holds a **simulated** dataset (fake genes on chr1–3) produced by `examples/make_example.py`. Causal variants were planted near TSSs, at splice sites, in open chromatin and on missense positions, so the report shows clear enrichments. It is for demos and testing only; do not draw biological conclusions from it.

### Web app reference annotations

Users upload only their credible sets. The reference annotations stay on the server and are listed in a JSON file, set via `FINEMAP_QC_REFS` (default: `examples/references.json`). Paths in the file are relative to the file itself:

```json
{"gtf": "gencode.v46.basic.annotation.gtf.gz",
 "bed": {"encode_ccre": "GRCh38-cCREs.bed"},
 "variant_set": {"missense": "missense.txt", "eqtlcat_rep": "eqtlcat_ge_rep.txt"},
 "score": {"spliceai": "spliceai.tsv"}}
```

```bash
FINEMAP_QC_REFS=/data/refs/references.json streamlit run app.py
```

Parsed GTF and BED files are cached in memory between runs, so only the first run pays the cost of loading the full GENCODE GTF (tens of seconds). Users can also upload extra annotation files in the sidebar.

## Input

A TSV with credible set ID, variant ID (`chr_pos_ref_alt`, GRCh38) and PIP. The tool picks the columns as follows:

* If there are columns named `cs_id`, `variant` and `pip`, it uses those. This means eQTL Catalogue `*.credible_sets.tsv.gz` files work as they are.
* Otherwise it uses the first three columns. A header line is optional.

Molecular trait: the tool uses `molecular_trait_id` if that column exists. Otherwise it takes `cs_id` and removes a trailing `_L<n>`.

Variant IDs in all inputs can use `_`, `:` or `-` as separators, with or without the `chr` prefix.

## Usage

```bash
finemap-qc --credible-sets QTD000001.credible_sets.tsv.gz \
  --gtf gencode.v46.basic.annotation.gtf.gz \
  --bed encode_ccre=GRCh38-cCREs.bed \
  --variant-set missense=missense_variants.txt \
  --variant-set finngen_rep=finngen_rep.txt \
  --variant-set eqtlcat_rep=eqtlcat_ge_rep.txt \
  --score spliceai=spliceai_max_delta.tsv \
  --out qc_out
```

(`python -m finemap_qc ...` also works without installing.)

| Option | Feature type | What it gives |
|---|---|---|
| `--gtf` | quantitative | `dist_tss` (distance to the nearest protein-coding transcript TSS) and `dist_splice_site` (distance to the nearest exon–intron boundary of a protein-coding transcript). Add `--canonical-only` to restrict to Ensembl canonical transcripts. |
| `--bed NAME=PATH` | binary | Whether the variant overlaps any interval, e.g. ENCODE cCREs or DNase peaks. |
| `--variant-set NAME=PATH` | binary | Whether the variant ID is in the list (first column), e.g. missense variants or replicating fine-mapped variants. |
| `--score NAME=PATH` | quantitative | Precomputed scores (variant, score). Variants missing from the table are excluded from that feature. |

### Variant groups

| Group | Definition |
|---|---|
| `high` | PIP > `--high-pip` (default 0.8) |
| `low_all` | PIP < `--low-pip` (default 0.05) |
| `low_anchored` | PIP < 0.05, in a credible set whose top variant has PIP > `--anchor-pip` (default 0.5) |

Each group counts unique variants. A variant that is high-PIP in any credible set is removed from both low groups.

### Output (`--out`)

* `report.html`: one self-contained page with all tables and plots.
* `descriptive_summary.tsv`:
  * number of credible sets and molecular traits
  * credible sets per trait
  * credible set size distribution
  * fraction of credible sets with max PIP > 0.5, > 0.8 and > 0.95
* `credible_sets.tsv`: trait, size and max PIP for each credible set.
* `binary_enrichment.tsv`: for each feature and comparison:
  * fraction of high and low variants that have the feature
  * fold enrichment with a 95% CI
  * Fisher odds ratio and p-value
* `quantitative_comparison.tsv`: medians, means, median difference and ratio, and Mann–Whitney p-value.
* `hist_<feature>.png`: overlaid histograms with lines at the medians. Distances are plotted as log10(d+1).
* `variant_annotations.tsv.gz`: all features for each variant, for your own downstream analysis.

## Reference data

### Downloaded automatically

```bash
python scripts/download_references.py --out-dir refs --exclude-study QTS000001
```
```bash
FINEMAP_QC_REFS=refs/references.json streamlit run app.py
```

The script fetches:

* **Gene model:** GENCODE basic GTF, release 46 by default (`--gencode-release`). About 50 MB.
* **Open chromatin:** ENCODE SCREEN cCREs, V3 registry (`--ccre-url` to use another file). The cCRE set is broad (about 8% of the genome), so tissue-matched DNase/ATAC peaks will discriminate better; pass those with `--bed` or add them to `references.json`.
* **Replicating eQTL Catalogue variants:** for each quantification method (`--eqtlcat-quant`, default `ge leafcutter`), the variants with PIP > 0.8 in at least 2 distinct studies. All datasets (tissues, conditions) of one study count as a single study. The script streams every credible set file from the EBI FTP and caches only the variants with PIP > 0.5 under `refs/eqtl_catalogue/`. The first run downloads all 280 `ge` and 113 `leafcutter` datasets, which takes a while. Re-runs with a different `--pip` (≥ 0.5), `--min-studies` or `--exclude-study` use the cache.
* **`references.json`:** the file the app and the CLI read. It also records which studies were excluded.

**Always pass `--exclude-study` with the eQTL Catalogue study ID of the dataset you are evaluating.** Otherwise the `eqtlcat_*` enrichments are circular. To evaluate a different study, re-run the script with that study's ID; it only takes seconds thanks to the cache.

The eQTL Catalogue has no MAJIQ quantifications, so those replicating sets can't be scripted.

### Manual steps

* **Missense:** run VEP once on the union of your input variants, keep the `missense_variant` ones, and add the list as a `variant_set`. This is cheaper and more exact than any genome-wide list.
* **Other fine-mapping resources (FinnGen, Open Targets, Tambets et al.):** download the files you want, then build a replicating set with `scripts/build_replicating_set.py`. Each input file counts as one study:
  ```bash
  python scripts/build_replicating_set.py --variant-col rsid --pip-col prob --out finngen_rep.txt *.SUSIE.snp.bgz
  ```
* **SpliceAI / ChromBPNet / AlphaGenome:** precompute these for the union of your input variants and add them as `score` tables (variant, score). The precomputed SpliceAI files require an Illumina login.

Add the files from the manual steps to `refs/references.json` (paths relative to that file). Re-running the download script keeps entries it didn't create.

## Caveats (read before interpreting numbers)

* **Circularity.** If the input is an eQTL Catalogue dataset, leave that dataset out of the eQTL Catalogue replicating set. Treat the same study and its tissues the same way, or the enrichment is inflated by construction.
* **P-values are anti-conservative.** Low-PIP variants within a credible set are in LD and spatially clustered, so they are not independent observations. Compare datasets by fold enrichment, not p-values. A proper CI would bootstrap over credible sets; that is not implemented.
* **`low_all` vs `low_anchored`.** `low_all` mixes in low-PIP variants from diffuse credible sets (weak or poorly resolved signals), which tend to sit in different genomic contexts. `low_anchored` is the fairer within-locus control, but it is much smaller.
* **Small high-PIP sets.** A small dataset may have only tens of high-PIP variants, which makes the CIs wide.
* **Allele matching.** Variant-set and score matching uses the exact `chrom_pos_ref_alt`. BED and distance features use position only. Indels are placed at their first base.
