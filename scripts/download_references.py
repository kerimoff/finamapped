#!/usr/bin/env python3
"""Download the reference annotations that can be fetched automatically, and write references.json.

    python scripts/download_references.py --out-dir refs --exclude-study QTS000001
    FINEMAP_QC_REFS=refs/references.json streamlit run app.py

What it fetches:
  * GENCODE basic GTF (GRCh38)                     -> dist_tss, dist_splice_site
  * ENCODE SCREEN cCREs (GRCh38)                   -> encode_ccre
  * eQTL Catalogue credible sets, per quant method -> eqtlcat_<quant>_rep
    (variants with PIP > --pip in >= --min-studies distinct studies; all
    datasets (tissues/conditions) of one study count as one study)

Not scripted (see README): missense (run VEP), SpliceAI and other scores,
FinnGen / Open Targets fine-mapping.

eQTL Catalogue files are streamed: only variants with PIP > 0.5 are cached per
dataset (refs/eqtl_catalogue/*.tsv), so re-running with another --pip >= 0.5,
--min-studies or --exclude-study is fast and needs no new downloads.

IMPORTANT: pass --exclude-study for the study you are evaluating, otherwise the
eqtlcat_* enrichments are circular.
"""
import argparse
import json
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from finemap_qc.io import normalize_variant_ids  # noqa: E402

GENCODE_URL = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
               "release_{r}/gencode.v{r}.basic.annotation.gtf.gz")
CCRE_URL = "https://downloads.wenglab.org/V3/GRCh38-cCREs.bed"
EQTLCAT_METADATA = ("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/"
                    "master/tabix/tabix_ftp_paths.tsv")
CACHE_PIP = 0.5


def retry(fn, tries=4):
    for i in range(tries):
        try:
            return fn()
        except Exception as e:  # network errors: back off and retry
            if i == tries - 1:
                raise
            print(f"    retrying after error: {e}")
            time.sleep(2 ** (i + 1))


def download(url, dest):
    if dest.exists():
        print(f"  exists, skipping: {dest}")
        return
    print(f"  downloading {url}")
    tmp = dest.with_name(dest.name + ".part")
    retry(lambda: urllib.request.urlretrieve(url, tmp))
    tmp.rename(dest)


def eqtl_catalogue(out_dir, metadata, quant_methods, pip, min_studies, exclude):
    """Return {quant_method: path of replicating variant list}."""
    if pip < CACHE_PIP:
        sys.exit(f"--pip must be >= {CACHE_PIP} (the cache threshold)")
    meta = retry(lambda: pd.read_csv(metadata, sep="\t", dtype=str))
    cache = out_dir / "eqtl_catalogue"
    cache.mkdir(exist_ok=True)
    out = {}
    for quant in quant_methods:
        rows = meta[(meta["quant_method"] == quant) & ~meta["study_id"].isin(exclude)]
        if rows.empty:
            print(f"  no eQTL Catalogue datasets for quant_method={quant!r}, skipping")
            continue
        print(f"  {quant}: {len(rows)} datasets from {rows['study_id'].nunique()} studies")
        study_hits = {}
        for i, r in enumerate(rows.itertuples(), 1):
            f = cache / f"{r.dataset_id}.tsv"
            if not f.exists():
                url = r.ftp_cs_path.replace("ftp://", "https://", 1)
                print(f"    [{i}/{len(rows)}] {r.dataset_id} ({r.study_label}, {r.sample_group})")
                df = retry(lambda: pd.read_csv(url, sep="\t", usecols=["variant", "pip"],
                                               compression="gzip" if url.endswith(".gz") else "infer"))
                df[df["pip"] > CACHE_PIP].to_csv(f, sep="\t", index=False)
            df = pd.read_csv(f, sep="\t")
            hits = set(normalize_variant_ids(df.loc[df["pip"] > pip, "variant"])["variant"].dropna())
            study_hits.setdefault(r.study_id, set()).update(hits)
        counts = pd.Series([v for hits in study_hits.values() for v in hits]).value_counts()
        keep = sorted(counts.index[counts >= min_studies])
        path = out_dir / f"eqtlcat_{quant}_rep.txt"
        path.write_text("".join(v + "\n" for v in keep))
        print(f"  {quant}: {len(keep):,} variants with PIP > {pip} in >= {min_studies} studies -> {path}")
        out[f"eqtlcat_{quant}_rep"] = path.name
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", default="refs")
    p.add_argument("--gencode-release", default="46")
    p.add_argument("--ccre-url", default=CCRE_URL)
    p.add_argument("--eqtlcat-metadata", default=EQTLCAT_METADATA, help="URL or path of tabix_ftp_paths.tsv")
    p.add_argument("--eqtlcat-quant", nargs="*", default=["ge", "leafcutter"],
                   help="quant methods to build replicating sets for (ge, exon, tx, txrev, leafcutter, ...)")
    p.add_argument("--exclude-study", nargs="*", default=[], metavar="QTS_ID",
                   help="eQTL Catalogue study IDs to leave out (the study you are evaluating!)")
    p.add_argument("--pip", type=float, default=0.8)
    p.add_argument("--min-studies", type=int, default=2)
    p.add_argument("--skip-gencode", action="store_true")
    p.add_argument("--skip-ccre", action="store_true")
    p.add_argument("--skip-eqtlcat", action="store_true")
    a = p.parse_args()

    out_dir = Path(a.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    refs_path = out_dir / "references.json"
    refs = json.loads(refs_path.read_text()) if refs_path.exists() else {}

    if not a.skip_gencode:
        print("GENCODE")
        gtf = out_dir / f"gencode.v{a.gencode_release}.basic.annotation.gtf.gz"
        download(GENCODE_URL.format(r=a.gencode_release), gtf)
        refs["gtf"] = gtf.name
    if not a.skip_ccre:
        print("ENCODE cCREs")
        bed = out_dir / Path(a.ccre_url).name
        download(a.ccre_url, bed)
        refs.setdefault("bed", {})["encode_ccre"] = bed.name
    if not a.skip_eqtlcat:
        print("eQTL Catalogue")
        sets = eqtl_catalogue(out_dir, a.eqtlcat_metadata, a.eqtlcat_quant, a.pip, a.min_studies,
                              set(a.exclude_study))
        refs.setdefault("variant_set", {}).update(sets)
        refs["eqtlcat_excluded_studies"] = sorted(a.exclude_study)

    refs_path.write_text(json.dumps(refs, indent=2) + "\n")
    print(f"\nWrote {refs_path}. Start the app with:\n  FINEMAP_QC_REFS={refs_path} streamlit run app.py")


if __name__ == "__main__":
    main()
