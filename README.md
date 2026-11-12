# finemap-qc

Quality metrics for fine-mapped credible sets in eQTL Catalogue format. The tool asks whether
high-PIP variants (PIP > 0.8) are enriched for features that true causal variants tend to have,
compared with low-PIP variants.

## Quick start

```bash
git clone https://github.com/kerimoff/finamapped.git && cd finamapped
git checkout hf-space
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

* If there are columns named `cs_id`, `variant` and `pip`, it uses those. Native eQTL Catalogue files are supported. Repeated rsID aliases require the explicit collapse option below.
* Headerless files use the first three columns. Named headers must contain `cs_id`, `variant` and `pip`.

PIPs must be finite and in [0, 1]; coordinates must be positive integers. Missing identifiers, conflicting trait/locus IDs within a credible set, and duplicate normalized credible-set/variant pairs are rejected. `--duplicate-policy collapse-identical` (also a checkbox in both web apps) permits identical analysis rows caused by multiple rsIDs; conflicting PIPs/traits/loci still fail. The report records removed rows. This does not validate the genome build or left-normalize indels.

Molecular trait: the tool uses `molecular_trait_id` if that column exists. Otherwise it takes `cs_id` and removes a trailing `_L<n>`.

Variant IDs in all inputs can use `_`, `:` or `-` as separators, with or without the `chr` prefix.

## Usage

```bash
finemap-qc --credible-sets QTD000001.credible_sets.tsv.gz \
  --duplicate-policy collapse-identical \
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
  * fold enrichment with a 95% cluster-bootstrap percentile CI
  * Fisher odds ratio and p-value
* `quantitative_comparison.tsv`: medians, means, median difference and ratio, 95% cluster-bootstrap median-difference CI, and exploratory Mann–Whitney p-value.
* `hist_<feature>.png`: overlaid histograms with lines at the medians. Distances are plotted as log10(d+1).
* `variant_clusters.tsv`: explicit variant-to-cluster assignments used for uncertainty.
* `variant_annotations.tsv.gz`: all features for each variant, for your own downstream analysis.

## Reference data

### Downloaded automatically

```bash
python scripts/download_references.py --out-dir refs
```
```bash
FINEMAP_QC_REFS=refs/references.json streamlit run app.py
```

The script fetches:

* **Gene model:** GENCODE basic GTF, release 46 by default (`--gencode-release`). About 50 MB.
* **Open chromatin:** ENCODE SCREEN cCREs, V3 registry (`--ccre-url` to use another file). The cCRE set is broad (about 8% of the genome), so tissue-matched DNase/ATAC peaks will discriminate better; pass those with `--bed` or add them to `references.json`.
* **Replicating eQTL Catalogue variants:** for each quantification method (`--eqtlcat-quant`, default `ge leafcutter`), the variants with PIP > 0.8 in at least 2 distinct studies. All datasets (tissues, conditions) of one study count as a single study. The script reads each credible set file from the EBI FTP and caches only the variants with PIP > 0.5 under `refs/eqtl_catalogue/`. The first run downloads all selected datasets in the metadata, which can take a while. Re-runs with a different `--pip` (≥ 0.5), `--min-studies` or `--exclude-study` use the cache.
* **`references.json`:** the file the app and the CLI read. It also records which studies were excluded.

finemap-qc is meant for fine-mapping results from outside the eQTL Catalogue (for example, a dataset submitted for inclusion), so the replicating sets use every study by default. To evaluate a dataset that is already in the eQTL Catalogue, re-run the script with `--exclude-study <QTS ID>` for that study; otherwise its `eqtlcat_*` enrichments are circular. Thanks to the cache this takes seconds.

The eQTL Catalogue has no MAJIQ quantifications, so those replicating sets can't be scripted.

### Manual steps

* **Missense:** run VEP once on the union of your input variants, keep the `missense_variant` ones, and add the list as a `variant_set`. This is cheaper and more exact than any genome-wide list.
* **Other fine-mapping resources (FinnGen, Open Targets, Tambets et al.):** download the files you want, then build a replicating set with `scripts/build_replicating_set.py`. Pass one study/cohort ID per file; repeated IDs are pooled (including multiple tissues or traits). Coordinate variant IDs are required, not rsIDs:
  ```bash
  python scripts/build_replicating_set.py cohortA.tsv cohortB.tsv --study-ids A B --out replicating.txt
  ```
* **SpliceAI / ChromBPNet / AlphaGenome:** precompute these for the union of your input variants and add them as `score` tables (variant, score). The precomputed SpliceAI files require an Illumina login.

Add the files from the manual steps to `refs/references.json` (paths relative to that file). Re-running the download script keeps entries it didn't create.

## Caveats (read before interpreting numbers)

* **Circularity (eQTL Catalogue inputs only).** External datasets are not affected. If the input is an eQTL Catalogue dataset, leave that dataset out of the eQTL Catalogue replicating set. Treat the same study and its tissues the same way, or the enrichment is inflated by construction.
* **P-values are anti-conservative.** Low-PIP variants within a credible set are in LD and spatially clustered, so they are not independent observations. Compare datasets by fold enrichment, not p-values. The CLI and web apps now bootstrap clusters for intervals; Fisher and Mann–Whitney p-values remain exploratory independence-based calculations, without multiple-testing correction.
* **`low_all` vs `low_anchored`.** `low_all` mixes in low-PIP variants from diffuse credible sets (weak or poorly resolved signals), which tend to sit in different genomic contexts. `low_anchored` restricts eligibility to anchored sets but is not a matched within-locus comparison; it can still be confounded.
* **Small high-PIP sets.** A small dataset may have only tens of high-PIP variants, which makes the CIs wide.
* **Allele matching.** Variant-set and score matching uses the exact `chrom_pos_ref_alt`. BED and distance features use position only. Indels are placed at their first base.

## Cluster uncertainty

The CLI and both web apps use 1,000 paired bootstrap draws (`--bootstrap-replicates`, minimum 100),
with `--bootstrap-seed 0`. By default, credible sets sharing a variant are joined into connected
components. An optional `locus_id` input column first joins credible sets assigned to the same locus;
shared variants then join any remaining overlapping components. Use curated independent locus/LD
blocks where possible. All variants in a component are resampled together, with replacement, retaining
draw multiplicity. Point estimates still count unique variants; duplicate row removal is separate.

Intervals describe the ratio of overlap fractions or the difference in medians. Each comparison/feature
uses its contributing clusters after missing annotations are excluded. Outputs include cluster counts,
seed, replicate count, valid draws and `ci_status`. Fewer than two contributing clusters in either group,
or any undefined bootstrap replicate (empty group/zero low overlap), suppresses the interval rather than
reporting a CI conditional on successful draws. Even two clusters are not enough for reliable inference;
inspect the counts. A missing interval is not evidence of no enrichment. Direct low-level `stats.py`
functions retain legacy independence-based tests/Wald intervals; the CLI replaces those intervals via
`uncertainty.cluster_intervals`.

Clusters based only on shared variants do not capture LD between disjoint sets. Bootstrap intervals
are not calibrated significance tests, do not correct confounding or annotation selection, and have
not undergone simulation-based coverage validation. Enrichment is not proof of causality or PIP calibration.

Reference preparation emits per-list `.provenance.json` files with input URLs, thresholds, study
exclusions and output checksums. Source URLs and cache version distinguish caches; a changed remote
file at the same URL still requires clearing its cache. Distinct study IDs can share participants, and
exact-variant recurrence across different traits is not necessarily replication of the same QTL.

See [real-data validation and limitations](docs/validation.md).
