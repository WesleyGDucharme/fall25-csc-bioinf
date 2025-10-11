#!/usr/bin/env bash
set -euo pipefail

# Run from this script's folder (week3/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# --- Activate venv if present (prefer ../.venv from repo root) ---
if [[ -f "../.venv/bin/activate" ]]; then
  source "../.venv/bin/activate"
elif [[ -f ".venv/bin/activate" ]]; then
  source ".venv/bin/activate"
fi

# Pick python executable (fallback to python3 if 'python' missing)
PYEXE="${PYEXE:-python}"
if ! command -v "$PYEXE" >/dev/null 2>&1; then
  PYEXE=python3
fi

# --- Python timing (Biotite baseline; no local Cython) ---
PY_OUT=""
if ! PY_OUT="$("$PYEXE" test/run_python_tests.py --just-ms --tests test/targeted_tests.py 2>&1)"; then
  echo "[ERROR] Python tests failed:"
  echo "$PY_OUT"
  exit 1
fi
PY_MS="$(echo "$PY_OUT" | tail -n1)"

# --- Codon timing ---
CODON_OUT=""
if ! CODON_OUT="$(
  export CODONPATH="$SCRIPT_DIR/code"
  cd "$SCRIPT_DIR/code" && codon run run_codon_tests.codon --just-ms 2>&1
)"; then
  echo "[ERROR] Codon tests failed:"
  echo "$CODON_OUT"
  exit 1
fi
CODON_MS="$(echo "$CODON_OUT" | tail -n1)"

# --- Required output ---
echo "Language    Runtime"
echo "-------------------"
printf "python      %sms\n" "$PY_MS"
printf "codon       %sms\n" "$CODON_MS"

