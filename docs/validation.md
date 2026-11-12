# Scientific correctness validation — 2026-11-12

Changes are based on `hf-space` commit `ad57fd9` (including the new Spaces frontend),
not the older default branch. The validation is a software and data-format check,
not a benchmark proving causal accuracy or calibrated PIPs.

## Public source data

- [eQTL Catalogue column definitions](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/tabix/Columns.md)
- [Catalogue metadata](https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv)
- [QTD000001 full credible sets](https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000001/QTD000001/QTD000001.credible_sets.tsv.gz): Alasoo_2018, naive macrophages, gene expression, QTS000001.
- [GENCODE v46 basic GRCh38 GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.basic.annotation.gtf.gz), with [GTF coordinate/attribute documentation](https://www.gencodegenes.org/pages/data_format.html).

SHA-256 checksums of downloaded inputs:

| Input | SHA-256 |
|---|---|
| QTD000001 credible sets | `f11924a944df7c32abfa928a3934ac25c7f0be9519d5b464d337935590f68021` |
| GENCODE v46 basic GTF | `d620d548ad23dad6c2d67486b7679d12f00f02f66dd5787183d8e6678dceb9b4` |
| Catalogue metadata | `20912adbeb700bb0224e73f1acba79921cb7afd4323955b00e4e9a73da49b7f4` |

`tests/data/QTD000001.one_cs.tsv` is the complete 56 source rows for
`ENSG00000230489_L1`, extracted from that public file. Its 53 distinct variants
match the source's `cs_size`; three extra rows carry alternative rsIDs. This small
fixture retains the native column order and is tested offline. No synthetic rows
are mixed into this fixture.

## Full-file run

Run with the two downloaded files in the current directory:

```bash
python -m finemap_qc --credible-sets QTD000001.credible_sets.tsv.gz \
  --duplicate-policy collapse-identical \
  --gtf gencode.v46.basic.annotation.gtf.gz --canonical-only \
  --bootstrap-replicates 1000 --bootstrap-seed 0 --out real-qc
```

Observed results:

- 71,634 source rows, 6,581 identical analysis duplicates removed, 65,053 retained.
- 1,331 credible sets, 1,213 molecular traits, 55,717 distinct variants.
- 99 high-PIP variants, 50,597 low-all variants, 1,078 low-anchored variants.
- 1,243 shared-variant connected components across the input.
- Canonical GENCODE selection yields 19,646 distinct chromosome/TSS positions and
  367,226 distinct chromosome/splice-boundary positions. These counts are not gene
  or transcript counts. Splice distance uses the last/first exonic base flanking
  an intron, in 1-based coordinates, on either strand; terminal transcript ends
  and single-exon transcripts are excluded.
- All four quantitative comparisons had 1,000/1,000 finite bootstrap draws.

| Feature | Comparison | Median high (bp) | Median low (bp) | Median difference 95% cluster interval (bp) |
|---|---|---:|---:|---|
| TSS distance | low_all | 9,496 | 19,167 | −13,836.225 to −3,797 |
| Splice distance | low_all | 1,860 | 4,533 | −4,370.5875 to −577.85 |
| TSS distance | low_anchored | 9,496 | 15,095.5 | −12,836.125 to 805.7 |
| Splice distance | low_anchored | 1,860 | 3,338.5 | −3,341.275 to 501.2 |

No biological significance claim is made from these exploratory results.

## Automated checks

Result: **40 tests passed**, including the optional Space test under Gradio 6.31.0.
Run `python -m pytest -q` (the Space test skips when Gradio is absent). Coverage includes invalid PIPs and identifiers, normalized
duplicates and conflicting trait IDs, native Catalogue columns, strict threshold
boundaries, canonical transcript-to-exon membership, both strands and single-exon
transcripts, shared-variant/locus components, bootstrap reproducibility and draw
multiplicity, undefined replicates, reference study pooling/exclusion, cache reuse,
provenance, local Space upload/error/report/ZIP endpoints, and CLI/report output including low-anchored group labels.

## Still unvalidated

- Frequentist coverage of bootstrap intervals, especially with few clusters; no
  simulation study or calibrated cluster-level p-values were added.
- True independence of clusters: disjoint credible sets can still be in LD.
  Supply curated `locus_id` values where available. Catalogue `region` intervals
  are not automatically interpreted as independent LD blocks.
- Independent-cohort, same-trait replication across the entire Catalogue. Tests
  verify study-ID pooling and exclusion, not participant overlap or biological
  concordance. This run did not download the full multi-study reference collection.
- ENCODE, FinnGen, Open Targets, Tambets, missense and prediction-score annotations
  were not validated with real files in this run. Genome-build matching, allele
  harmonization and indel left-normalization remain the user's responsibility.
- Live Spaces deployment, remote reference-dataset availability, memory/queue load
  and hosted authentication are separate from local CLI/runtime checks.
