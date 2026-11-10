"""finemap-qc web app for Hugging Face Spaces (Gradio SDK, custom frontend via gr.Server).

Reference annotations are downloaded at start-up from the dataset repo named in
$FINEMAP_QC_REFS_REPO (it must contain a references.json).  Without it, only the
simulated example annotations are available.

Runs execute in a background thread; the page polls /api/runs/<id> for status.
"""
import io
import json
import os
import re
import shutil
import tempfile
import threading
import time
import traceback
import uuid
import zipfile
from pathlib import Path

try:  # ZeroGPU needs at least one decorated function; this one is never called (CPU-only app)
    import spaces

    @spaces.GPU(duration=1)
    def _noop():
        pass
except ImportError:  # running locally
    pass

from fastapi import File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from gradio import Server

from finemap_qc import annotate, cli, io as fio

ROOT = Path(__file__).resolve().parent
EXAMPLE_REFS = ROOT / "examples" / "references.json"
EXAMPLE_CS = ROOT / "examples" / "credible_sets.tsv.gz"
REFS_REPO = os.environ.get("FINEMAP_QC_REFS_REPO", "").strip()
WORK = Path(tempfile.gettempdir()) / "finemap_qc_runs"
WORK.mkdir(exist_ok=True)
MAX_AGE = 2 * 3600  # seconds a finished run is kept

state = {"refs": None, "refs_status": "none" if not REFS_REPO else "loading", "refs_error": None}
runs = {}
run_lock = threading.Semaphore(1)  # one analysis at a time keeps memory predictable


def _ref_args(path):
    cfg = json.loads(Path(path).read_text())
    res = lambda p: str((Path(path).parent / p).resolve())  # noqa: E731
    args = ["--gtf", res(cfg["gtf"])] if cfg.get("gtf") else []
    feats = []
    if cfg.get("gtf"):
        feats += [{"name": "dist_tss", "kind": "quantitative", "source": cfg["gtf"]},
                  {"name": "dist_splice_site", "kind": "quantitative", "source": cfg["gtf"]}]
    for kind in ("bed", "variant_set", "score"):
        for name, p in cfg.get(kind, {}).items():
            args += [f"--{kind.replace('_', '-')}", f"{name}={res(p)}"]
            feats.append({"name": name, "kind": "quantitative" if kind == "score" else "binary", "source": p})
    return {"args": args, "features": feats, "excluded_studies": cfg.get("eqtlcat_excluded_studies", []),
            "gtf": res(cfg["gtf"]) if cfg.get("gtf") else None,
            "beds": [res(p) for p in cfg.get("bed", {}).values()]}


def _warm(refs):
    """Parse the big reference files once so the first run is fast (both functions are lru_cached)."""
    if refs["gtf"]:
        annotate.load_gene_features(refs["gtf"], False)
    for b in refs["beds"]:
        fio.read_bed(b)


def _load_refs():
    try:
        from huggingface_hub import snapshot_download
        local = snapshot_download(REFS_REPO, repo_type="dataset", local_dir=str(ROOT / "refs"),
                                  allow_patterns=["*.json", "*.gz", "*.bed", "*.txt", "*.tsv"])
        state["refs"] = _ref_args(Path(local) / "references.json")
        state["refs_status"] = "warming"
        _warm(state["refs"])
        state["refs_status"] = "ready"
    except Exception as e:  # surface in the UI instead of crashing the Space
        state["refs_status"], state["refs_error"] = "error", f"{type(e).__name__}: {e}"
        traceback.print_exc()


if REFS_REPO:
    threading.Thread(target=_load_refs, daemon=True).start()
EXAMPLE = _ref_args(EXAMPLE_REFS)


def _cleanup():
    now = time.time()
    for rid, r in list(runs.items()):
        if r.get("finished") and now - r["finished"] > MAX_AGE:
            shutil.rmtree(r["dir"], ignore_errors=True)
            runs.pop(rid, None)


def _safe(name, default):
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(name or default).name).strip("._") or default
    return stem[:80]


def _save(upload, folder, stem):
    suffix = "".join(Path(upload.filename or "").suffixes[-2:]) or ".tsv"
    folder.mkdir(exist_ok=True)
    path = folder / f"{stem}{suffix}"
    with open(path, "wb") as fh:
        shutil.copyfileobj(upload.file, fh)
    return path


def _execute(rid, argv):
    r = runs[rid]
    with run_lock:
        r["status"], r["started"] = "running", time.time()
        try:
            cli.main(argv)
            out = Path(r["dir"]) / "out"
            buf = io.BytesIO()
            with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
                for f in sorted(out.iterdir()):
                    if f.suffix != ".svg":
                        z.write(f, f.name)
            (Path(r["dir"]) / "results.zip").write_bytes(buf.getvalue())
            r["status"] = "done"
        except SystemExit as e:  # cli uses sys.exit for input errors
            r["status"], r["error"] = "error", str(e)
        except Exception as e:
            r["status"], r["error"] = "error", (
                f"Could not process the input ({type(e).__name__}: {e}). Check that the credible sets file "
                "is tab-separated with cs_id, variant and pip columns, and that annotation files match "
                "the formats described above.")
            traceback.print_exc()
        finally:
            r["finished"] = time.time()


app = Server(title="finemap-qc")


@app.get("/", response_class=HTMLResponse)
async def homepage():
    return (ROOT / "index.html").read_text()


@app.get("/api/config")
def config():
    refs = state["refs"]
    return {"refs_status": state["refs_status"], "refs_error": state["refs_error"],
            "refs_repo": REFS_REPO or None,
            "features": refs["features"] if refs else [],
            "excluded_studies": refs["excluded_studies"] if refs else [],
            "example_features": EXAMPLE["features"], "busy": run_lock._value == 0}


@app.post("/api/run")
def start_run(credible_sets: UploadFile = File(None), use_example: bool = Form(False),
              high_pip: float = Form(0.8), low_pip: float = Form(0.05), anchor_pip: float = Form(0.5),
              canonical_only: bool = Form(False), use_refs: bool = Form(True),
              bed: list[UploadFile] = File(None), variant_set: list[UploadFile] = File(None),
              score: list[UploadFile] = File(None)):
    _cleanup()
    if not 0 <= low_pip < high_pip <= 1 or not 0 <= anchor_pip <= 1:
        raise HTTPException(400, "Thresholds must satisfy 0 ≤ low PIP < high PIP ≤ 1.")
    if not use_example and credible_sets is None:
        raise HTTPException(400, "Upload a credible sets file or run the example.")
    rid = uuid.uuid4().hex[:12]
    d = WORK / rid
    d.mkdir()
    if use_example:
        cs, name, ref = EXAMPLE_CS, EXAMPLE_CS.name, EXAMPLE["args"]
    else:
        name = _safe(credible_sets.filename, "credible_sets.tsv")
        cs = _save(credible_sets, d, name.split(".")[0])
        if use_refs and state["refs_status"] in ("loading", "warming"):
            raise HTTPException(409, "Reference annotations are still loading. Try again in a minute.")
        ref = state["refs"]["args"] if use_refs and state["refs"] else []
    argv = ["--credible-sets", str(cs), "--out", str(d / "out"), "--high-pip", str(high_pip),
            "--low-pip", str(low_pip), "--anchor-pip", str(anchor_pip)] + list(ref)
    if canonical_only:
        argv.append("--canonical-only")
    if not use_example:
        for flag, files in (("--bed", bed), ("--variant-set", variant_set), ("--score", score)):
            for f in files or []:
                if not f.filename:
                    continue
                fname = re.sub(r"\W+", "_", Path(f.filename).name.split(".")[0]).strip("_") or "feature"
                argv += [flag, f"{fname}={_save(f, d / flag.strip('-'), fname)}"]
    queued = run_lock._value == 0
    runs[rid] = {"dir": str(d), "status": "queued" if queued else "starting", "name": name,
                 "created": time.time(), "error": None, "finished": None}
    threading.Thread(target=_execute, args=(rid, argv), daemon=True).start()
    return {"id": rid}


def _run(rid):
    r = runs.get(rid)
    if not r:
        raise HTTPException(404, "Run not found (results are kept for 2 hours).")
    return r


@app.get("/api/runs/{rid}")
def run_status(rid: str):
    r = _run(rid)
    t0 = r.get("started") or r["created"]
    return {"status": r["status"], "error": r["error"], "name": r["name"],
            "elapsed": round((r["finished"] or time.time()) - t0, 1)}


@app.get("/runs/{rid}/report.html")
def report(rid: str):
    p = Path(_run(rid)["dir"]) / "out" / "report.html"
    if not p.exists():
        raise HTTPException(404, "Report not ready.")
    return HTMLResponse(p.read_text())


@app.get("/runs/{rid}/results.zip")
def results(rid: str):
    p = Path(_run(rid)["dir"]) / "results.zip"
    if not p.exists():
        raise HTTPException(404, "Results not ready.")
    return FileResponse(p, media_type="application/zip", filename="finemap_qc_results.zip")


@app.get("/runs/{rid}/report-download")
def report_download(rid: str):
    p = Path(_run(rid)["dir"]) / "out" / "report.html"
    if not p.exists():
        raise HTTPException(404, "Report not ready.")
    return Response(p.read_bytes(), media_type="text/html",
                    headers={"Content-Disposition": 'attachment; filename="finemap_qc_report.html"'})


demo = app  # the HF runtime looks for `demo`
if __name__ == "__main__":
    demo.launch(ssr_mode=False)
