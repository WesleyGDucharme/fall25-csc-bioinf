set -euo pipefail

# ---------- build codon test once ----------
BIN_DIR=".bin"
BIN="$BIN_DIR/codon_test"
mkdir -p "$BIN_DIR"
# Rebuild if sources are newer than the binary
if [[ ! -x "$BIN" || code/codon_test.codon -nt "$BIN" || code/codon_align.codon -nt "$BIN" ]]; then
  echo "[build] codon -> $BIN"
  codon build -release -o "$BIN" code/codon_test.codon
fi

# ---------- runners ----------
PY="python3 test/py_test.py"
CODON="$BIN"   # use the compiled binary

# For MT pair, affine can OOM; keep a safe default. Override when needed:
MT_AFFINE_FLAGS="${MT_AFFINE_FLAGS:---limit 3000}"
MT_LINEAR_FLAGS="${MT_LINEAR_FLAGS:-}"

# method sets (and nice labels)
methods=(global-linear local-linear fitting-affine global-affine)
label_of() {
  case "$1" in
    global-linear)  echo "globalL" ;;
    local-linear)   echo "local" ;;
    fitting-affine) echo "fittingA" ;;
    global-affine)  echo "globalA" ;;
    *)              echo "$1" ;;
  esac
}

# time helper -> integer ms
run_ms() {
  local cmd="$*"
  local sec
  { TIMEFORMAT=%R; sec=$( { time bash -c "$cmd" >/dev/null; } 2>&1 ); } 2>/dev/null
  awk -v s="$sec" 'BEGIN{printf("%.0f", s*1000)}'
}

print_header() {
  printf "%-18s %-10s %s\n" "Method" "Language" "Runtime"
  printf -- "--------------------------------------\n"
}

run_case() {
  local label="$1" lang="$2" method="$3" extra="$4"
  local ms
  if [[ "$lang" == "python" ]]; then
    ms=$(run_ms "$PY --method $method $extra")
  else
    ms=$(run_ms "$CODON --method $method $extra")
  fi
  printf "%-18s %-10s %sms\n" "$label" "$lang" "$ms"
}

print_header

# ---------- MT human vs orang (human->orang) ----------
for lang in python codon; do
  run_case "globalL-mt_human" "$lang" "global-linear"  "--query data/MT-human.fa --target data/MT-orang.fa --fasta $MT_LINEAR_FLAGS"
  run_case "local-mt_human"   "$lang" "local-linear"   "--query data/MT-human.fa --target data/MT-orang.fa --fasta $MT_LINEAR_FLAGS"
  run_case "fittingA-mt_human" "$lang" "fitting-affine" "--query data/MT-human.fa --target data/MT-orang.fa --fasta $MT_AFFINE_FLAGS"
  run_case "globalA-mt_human"  "$lang" "global-affine"  "--query data/MT-human.fa --target data/MT-orang.fa --fasta $MT_AFFINE_FLAGS"
done

# ---------- q1..q5 vs t1..t5 ----------
for i in 1 2 3 4 5; do
  for lang in python codon; do
    for m in "${methods[@]}"; do
      lbl="$(label_of "$m")-q$i"
      run_case "$lbl" "$lang" "$m" "--query data/q1.fa --target data/t1.fa --fasta --index $i"
    done
  done
done
