#!/usr/bin/env bash
set -euo pipefail

export CODON_PYTHON=/usr/lib/x86_64-linux-gnu/libpython3.12.so.1.0

# Directories (adjust if your layout differs)
CODE_DIR="$(cd "$(dirname "$0")"/code && pwd)"    # ../code from week1/
DATA_DIR="$(cd "$(dirname "$0")"/data && pwd)"    # ../data from week1/

# Datasets to evaluate
DATASETS=(data1 data2 data3 data4)

# For data4 we need a big stack (both Python & Codon). Others are fine.
NEEDS_BIG_STACK=(data4)

# Format seconds -> H:MM:SS
format_time () {
  local secs="$1"
  local h=$((secs/3600))
  local m=$(((secs%3600)/60))
  local s=$((secs%60))
  printf "%d:%02d:%02d" "$h" "$m" "$s"
}

# Compute N50 from a FASTA file that contains contigs
n50_from_fasta () {
  python3 - "$1" <<'PY'
import sys
path = sys.argv[1]
lengths = []
L = 0
with open(path, "r") as f:
    cur = []
    for line in f:
        if line.startswith(">"):
            if cur:
                seq = "".join(cur)
                lengths.append(len(seq))
                L += len(seq)
                cur = []
        else:
            cur.append(line.strip())
    if cur:
        seq = "".join(cur)
        lengths.append(len(seq))
        L += len(seq)
if not lengths:
    print(0)
    sys.exit(0)
lengths.sort(reverse=True)
half = L/2
cum = 0
for l in lengths:
    cum += l
    if cum >= half:
        print(l)
        break
PY
}

# Header
printf "Dataset\tLanguage\tRuntime\t\tN50\n"
printf -- "-------------------------------------------------------------------------------------------------------\n"

for ds in "${DATASETS[@]}"; do
  for lang in python codon; do
    out_fasta="${DATA_DIR}/${ds}/contig.fasta"
    # Clean previous contig (so we don't accidentally read stale results)
    rm -f "$out_fasta"

    start=$(date +%s)

    # Stack boost for data4 (Codon and Python)
    if printf "%s\n" "${NEEDS_BIG_STACK[@]}" | grep -qx "$ds"; then
      ulimit -s 8192000
    fi

    if [[ "$lang" == "python" ]]; then
      # CPython run
      ( cd "$CODE_DIR" && python3 main.py "${DATA_DIR}/${ds}" >/dev/null )
    else
      # Codon run
      ( cd "$CODE_DIR" && codon run -release main.py "${DATA_DIR}/${ds}" >/dev/null )
    fi

    end=$(date +%s)
    dur=$(( end - start ))
    pretty=$(format_time "$dur")

    # Compute N50 from the freshly written contig.fasta
    if [[ -f "$out_fasta" ]]; then
      n50=$(n50_from_fasta "$out_fasta")
    else
      n50=0
    fi

    printf "%s\t%s\t\t%s\t\t%s\n" "$ds" "$lang" "$pretty" "$n50"
  done
done
