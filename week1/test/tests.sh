#!/usr/bin/env bash
set -euo pipefail

# Resolve paths relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE_DIR="${CODE_DIR:-"$SCRIPT_DIR/../code"}"
DATA_DIR="${DATA_DIR:-"$SCRIPT_DIR/../data"}"

# Make diffs predictable
export LC_ALL=C

# Try to auto-set CODON_PYTHON if not provided
if [[ -z "${CODON_PYTHON:-}" ]]; then
  CAND=$(python3 - <<'PY'
import sysconfig, glob, os
libdir = sysconfig.get_config_var("LIBDIR") or "/usr/lib/x86_64-linux-gnu"
cands = sorted(glob.glob(os.path.join(libdir, "libpython*.so*")))
print(cands[-1] if cands else "")
PY
)
  if [[ -n "$CAND" ]]; then
    export CODON_PYTHON="$CAND"
  fi
fi

printf "Python: %s\n" "$(python3 -V 2>&1)"
printf "Codon : %s\n" "$(codon --version 2>&1 || true)"
printf "CODON_PYTHON: %s\n" "${CODON_PYTHON:-"(not set)"}"
echo

# quick sanity check for Codon->Python interop (non-fatal)
cat >"$SCRIPT_DIR/.codon_sanity.codon" <<'CODON'
from python import sys
print("Codon sees Python:", sys.version.split()[0])
CODON
if ! codon run "$SCRIPT_DIR/.codon_sanity.codon" >/dev/null 2>&1; then
  echo "⚠️  Codon Python interop sanity check failed. If Codon runs differ, set CODON_PYTHON explicitly."
fi
rm -f "$SCRIPT_DIR/.codon_sanity.codon"

# Clean old results
rm -f "$SCRIPT_DIR"/py_data*.out "$SCRIPT_DIR"/co_data*.out

# Helper to maybe lift stack for data4
maybe_stack() {
  local ds="$1"
  if [[ "$ds" == "data4" ]]; then
    ulimit -s 8192000 || true
  fi
}

# Run one pair (Python & Codon) with timing
run_pair() {
  local ds="$1"
  echo "=== $ds ==="
  maybe_stack "$ds"

  /usr/bin/time -f "PY time: %E" \
    python3 "$CODE_DIR/main.py" "$DATA_DIR/$ds" | tee "$SCRIPT_DIR/py_${ds}.out" >/dev/null

  maybe_stack "$ds"
  /usr/bin/time -f "CO time: %E" \
    codon run -release "$CODE_DIR/main.py" "$DATA_DIR/$ds" | tee "$SCRIPT_DIR/co_${ds}.out" >/dev/null

  echo
}

echo "Running tests (timed)…"
run_pair data1
run_pair data2
run_pair data3
run_pair data4

echo "Comparing outputs (Python vs Codon)"
for d in data1 data2 data3 data4; do
  echo "--- $d ---"
  if diff -u "$SCRIPT_DIR/py_${d}.out" "$SCRIPT_DIR/co_${d}.out"; then
    echo "$d: ✅ match"
  else
    echo "$d: ❌ mismatch"
    # show a quick focused view of first difference in lengths (optional)
    echo "First differing lines (if any):"
    diff -u "$SCRIPT_DIR/py_${d}.out" "$SCRIPT_DIR/co_${d}.out" | sed -n '1,40p'
  fi
  echo
done

