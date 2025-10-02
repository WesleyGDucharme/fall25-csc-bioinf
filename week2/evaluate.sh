set -euo pipefail

# Where am I?
ROOT="$(cd "$(dirname "$0")" && pwd)"
CODE="$ROOT/code"
TEST="$ROOT/test"
DATA="$ROOT/data"

# Optional: ensure Codon is on PATH (only if user installed in ~/.codon)
export PATH="$HOME/.codon/bin:$PATH"

echo "==> Python venv: $(python -V || true)"
echo "==> Codon: $(which codon || echo 'codon not found')"

# Build Codon worker (library backend for TRVIZ_IMPL=codon)
echo
echo "==> Building Codon worker…"
codon build -release -o "$CODE/trviz_codon_worker" "$CODE/trviz_codon/worker.codon"

#  Build Codon CLI (standalone pipeline)
echo
echo "==> Building Codon CLI…"
codon build -release -o "$CODE/trviz_codon_cli" "$CODE/trviz_codon/main.codon"

# Run Python (pure) tests
echo
echo "==> Running tests (Python backend)…"
( cd "$TEST" && PYTHONPATH=../code TRVIZ_IMPL=python pytest -q test_decomposer.py test_refinement.py )

# Run tests with Codon backend (Python -> worker bridge)
echo
echo "==> Running tests (Codon backend via worker)…"
( cd "$TEST" && PYTHONPATH=../code TRVIZ_IMPL=codon pytest -q test_decomposer.py test_refinement.py )

# Codon port smoke tests unit
echo
echo "==> Codon port smoke tests…"
export PYTHONPATH="$ROOT/code"
export TRVIZ_IMPL=codon
pytest -q "$ROOT/test/test_codon_port_smoke.py"

# CLI smoke test with Codon main
echo
echo "==> CLI smoke test (Codon main)…"

# default OUT_TSV if not provided by env
OUT_TSV="${OUT_TSV:-$ROOT/trviz_decomposed.tsv}"

# Run the Codon CLI to generate output
"$CODE/trviz_codon_cli" \
  --input "$DATA/demo.fasta" \
  --motifs "ACT,TT,G" \
  --aligner none \
  --out "$OUT_TSV"

echo "==> CLI wrote: $OUT_TSV"
head -n 5 "$OUT_TSV" || true

echo
echo "==> All done. ✅"
