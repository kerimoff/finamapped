#!/usr/bin/env bash
# Deploy the web app to a Hugging Face Space.
#   scripts/deploy_space.sh <user>/<space> [<user>/<refs-dataset>]
# With a dataset ID, refs/ (built by scripts/download_references.py) is uploaded to it first,
# and the Space is pointed at it. Requires `hf auth login`.
set -euo pipefail
SPACE=${1:?usage: deploy_space.sh <user>/<space> [<user>/<refs-dataset>]}
REFS=${2:-}
cd "$(dirname "$0")/.."

if [ -n "$REFS" ]; then
  [ -f refs/references.json ] || { echo "refs/references.json missing: run scripts/download_references.py first" >&2; exit 1; }
  hf repos create "$REFS" --type dataset --exist-ok
  hf upload "$REFS" refs . --repo-type dataset --exclude "eqtl_catalogue/**" --commit-message "Update reference annotations"
fi

hf repos create "$SPACE" --type space --space-sdk gradio --flavor zero-a10g --exist-ok
STAGE=$(mktemp -d)
trap 'rm -rf "$STAGE"' EXIT
cp space/app.py space/index.html space/requirements.txt space/README.md "$STAGE"/
cp -R finemap_qc "$STAGE"/finemap_qc
mkdir "$STAGE"/examples
cp examples/credible_sets.tsv.gz examples/genes.gtf.gz examples/open_chromatin.bed examples/missense.txt \
   examples/spliceai.tsv examples/references.json "$STAGE"/examples/
hf spaces variables add "$SPACE" --env GRADIO_SSR_MODE=false
[ -n "$REFS" ] && hf spaces variables add "$SPACE" --env FINEMAP_QC_REFS_REPO="$REFS"
hf upload "$SPACE" "$STAGE" . --repo-type space --exclude "**/__pycache__/**" --commit-message "Deploy finemap-qc"
echo "Deployed: https://huggingface.co/spaces/$SPACE"
