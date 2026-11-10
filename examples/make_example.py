#!/usr/bin/env python3
"""Generate the synthetic example dataset in this folder.

Everything here is SIMULATED (fake genes on chr1-3), purely to demo the tool.
Causal variants are placed preferentially near TSSs, at splice sites, in open
chromatin and on missense positions, so the report shows real enrichments:

    credible_sets.tsv.gz   eQTL Catalogue column layout (molecular_trait_id, cs_id, variant, pip, ...)
    genes.gtf.gz           protein-coding transcripts/exons
    open_chromatin.bed     peaks (TSS-proximal + random)
    missense.txt           variant IDs
    spliceai.tsv           variant, score
    references.json        annotation set used by the web app / CLI
"""
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
rng = np.random.default_rng(42)
N_GENES_PER_CHROM, GENE_SPACING = 150, 300_000

gtf, peaks, cs_rows, missense, scores = [], [], [], set(), {}
attrs = 'gene_id "{g}"; transcript_id "{t}"; gene_type "protein_coding"; transcript_type "protein_coding";'

for chrom in ("1", "2", "3"):
    for i in range(N_GENES_PER_CHROM):
        gene = f"ENSG{chrom}{i:05d}"
        strand = rng.choice(["+", "-"])
        start = 1_000_000 + i * GENE_SPACING
        n_exons = int(rng.integers(3, 12))
        exon_starts = np.sort(rng.choice(np.arange(start, start + 60_000, 200), n_exons, replace=False))
        exons = [(s, s + int(rng.integers(80, 250))) for s in exon_starts]
        end = exons[-1][1]
        a = attrs.format(g=gene, t=gene.replace("ENSG", "ENST"))
        gtf.append(f"chr{chrom}\tSIM\ttranscript\t{start}\t{end}\t.\t{strand}\t.\t{a} tag \"Ensembl_canonical\";")
        gtf += [f"chr{chrom}\tSIM\texon\t{s}\t{e}\t.\t{strand}\t.\t{a}" for s, e in exons]
        tss = start if strand == "+" else end
        peaks.append((chrom, tss - 1000, tss + 1000))
        peaks += [(chrom, p, p + 500) for p in rng.integers(start - 100_000, end + 100_000, 15)]
        introns = [e for _, e in exons[:-1]] + [s for s, _ in exons[1:]]

        for L in range(1, 1 + int(rng.choice([1, 1, 1, 2, 3]))):
            r = rng.random()  # where the causal variant sits
            if r < 0.35:
                causal = tss + int(rng.normal(0, 1500))
            elif r < 0.5:
                causal = int(rng.choice(introns)) + int(rng.integers(-5, 6))
            elif r < 0.6:
                s, e = exons[int(rng.integers(n_exons))]
                causal = int(rng.integers(s, e))
            else:
                causal = tss + int(rng.integers(-150_000, 150_000))
            n = int(rng.geometric(1 / 25))
            others = causal + rng.normal(0, 25_000, n - 1).astype(int)
            positions = np.concatenate([[causal], others])
            # resolved signals put most of the PIP on the causal variant; diffuse ones spread it
            pips = rng.dirichlet(np.ones(n) * (0.05 if rng.random() < 0.5 else 1.0))
            pips = np.sort(pips)[::-1]
            for j, (pos, pip) in enumerate(zip(positions, pips)):
                ref, alt = rng.choice(list("ACGT"), 2, replace=False)
                vid = f"chr{chrom}_{pos}_{ref}_{alt}"
                cs_rows.append((gene, gene, f"{gene}_L{L}", vid, f"rs{rng.integers(1e8)}", n, pip))
                at_splice = min(abs(pos - x) for x in introns) <= 5
                in_exon = any(s <= pos <= e for s, e in exons)
                if in_exon and (j == 0 or rng.random() < 0.3):
                    missense.add(vid)
                scores[vid] = round(float(rng.beta(5, 3) if at_splice and j == 0 else rng.beta(1, 40)), 3)

cs = pd.DataFrame(cs_rows, columns=["molecular_trait_id", "gene_id", "cs_id", "variant", "rsid", "cs_size", "pip"])
cs["pip"] = cs["pip"].round(5)
cs.to_csv(HERE / "credible_sets.tsv.gz", sep="\t", index=False)
with gzip.open(HERE / "genes.gtf.gz", "wt") as fh:
    fh.write("\n".join(gtf) + "\n")
pd.DataFrame(peaks).to_csv(HERE / "open_chromatin.bed", sep="\t", header=False, index=False)
(HERE / "missense.txt").write_text("".join(v + "\n" for v in sorted(missense)))
pd.Series(scores, name="spliceai").rename_axis("variant").to_csv(HERE / "spliceai.tsv", sep="\t")
(HERE / "references.json").write_text(json.dumps({
    "gtf": "genes.gtf.gz",
    "bed": {"open_chromatin": "open_chromatin.bed"},
    "variant_set": {"missense": "missense.txt"},
    "score": {"spliceai": "spliceai.tsv"},
}, indent=2) + "\n")
print(f"{cs['cs_id'].nunique()} credible sets, {len(cs)} rows written to {HERE}")
