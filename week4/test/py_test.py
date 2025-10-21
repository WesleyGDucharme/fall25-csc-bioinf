# Python version of tests

"""
Python alignment test harness.

Examples:
  # q1 vs t1 (first sequences)
  python3 test/py_test.py --method global-linear --query data/q1.fa --target data/t1.fa --fasta

  # q3 vs t3 by index (1-based)
  python3 test/py_test.py --method local-linear --query data/q1.fa --target data/t1.fa --fasta --index 3

  # q5 vs t5 by header token
  python3 test/py_test.py --method global-affine --query data/q1.fa --target data/t1.fa --fasta --qid q5 --tid t5
"""
import argparse
import sys
from pathlib import Path

# Make ../code importable
REPO_ROOT = Path(__file__).resolve().parents[1]
CODE_DIR = REPO_ROOT / "code"
sys.path.insert(0, str(CODE_DIR))

try:
    import py_align  # expected to define the four alignment functions
except Exception as e:
    print(f"IMPORT_ERROR\t{e}", file=sys.stderr)
    sys.exit(2)


def read_all_fasta(path: Path):
    """Return list of (header, seq) from a FASTA file. Uppercases the sequence."""
    records = []
    header = None
    seq_parts = []
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts).upper()))
                header = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_parts).upper()))
    return records


def select_seq(records, index=None, wanted_id=None, label="query"):
    """
    Pick a sequence from FASTA records:
      - If wanted_id (header token) is given, use it (exact match on first token).
      - Else if index (1-based) is given, pick that record.
      - Else pick the first record.
    """
    if wanted_id:
        for h, s in records:
            if h.split()[0] == wanted_id:
                return s
        print(f"ARG_ERROR\t{label} id '{wanted_id}' not found", file=sys.stderr)
        sys.exit(2)
    if index is not None:
        i = index - 1
        if i < 0 or i >= len(records):
            print(f"ARG_ERROR\t{label} index {index} out of range (n={len(records)})", file=sys.stderr)
            sys.exit(2)
        return records[i][1]
    if not records:
        print(f"ARG_ERROR\t{label} FASTA is empty", file=sys.stderr)
        sys.exit(2)
    return records[0][1]


def read_raw_or_fasta(path: Path, is_fasta: bool, index=None, wanted_id=None, label="query") -> str:
    if is_fasta:
        recs = read_all_fasta(path)
        return select_seq(recs, index=index, wanted_id=wanted_id, label=label)
    return path.read_text(encoding="utf-8").strip().upper()

def slice_1based(seq: str, spec: str) -> str:
    if not spec:
        return seq
    try:
        a, b = spec.split(":")
        i = int(a) if a else 1
        j = int(b) if b else len(seq)
    except Exception:
        print(f"ARG_ERROR\tbad slice '{spec}'", file=sys.stderr); sys.exit(2)
    i = max(1, i); j = min(len(seq), j)
    if j < i: return ""
    return seq[i-1:j]

def main() -> int:
    p = argparse.ArgumentParser(description="Python alignment test harness")
    p.add_argument("--method", required=True,
                   choices=["global-linear", "local-linear", "fitting-affine", "global-affine"])
    p.add_argument("--query", required=True, help="Path to query (FASTA or raw)")
    p.add_argument("--target", required=True, help="Path to target (FASTA or raw)")
    p.add_argument("--fasta", action="store_true", help="Interpret inputs as FASTA")
    p.add_argument("--index", type=int, help="1-based pair index from each FASTA (qN vs tN)")
    p.add_argument("--qid", type=str, help="query FASTA header token to select (e.g., q3)")
    p.add_argument("--tid", type=str, help="target FASTA header token to select (e.g., t3)")
    p.add_argument("--limit", type=int, help="truncate both sequences to first N bases")
    p.add_argument("--slice-q", type=str, help="1-based inclusive slice for query, e.g. 1:4000")
    p.add_argument("--slice-t", type=str, help="1-based inclusive slice for target, e.g. 1:4000")
    args = p.parse_args()

    # Precedence: qid/tid > index > first sequence
    q = read_raw_or_fasta(Path(args.query), args.fasta,
                          index=(None if args.qid else args.index),
                          wanted_id=args.qid, label="query")
    t = read_raw_or_fasta(Path(args.target), args.fasta,
                          index=(None if args.tid else args.index),
                          wanted_id=args.tid, label="target")

    # apply slices/limit (in this order: slice then limit)
    q = slice_1based(q, args.slice_q or "")
    t = slice_1based(t, args.slice_t or "")
    if args.limit:
        q = q[:args.limit]
        t = t[:args.limit]

    method_map = {
        "global-linear": "global_align_linear",
        "local-linear": "local_align_linear",
        "fitting-affine": "fitting_align_affine",
        "global-affine": "global_align_affine",
    }
    func_name = method_map[args.method]
    if not hasattr(py_align, func_name):
        print(f"NOT_IMPLEMENTED\t{func_name}", file=sys.stderr)
        return 2

    func = getattr(py_align, func_name)
    try:
        score, aq, at = func(q, t)
    except Exception as e:
        print(f"RUNTIME_ERROR\t{e}", file=sys.stderr)
        return 2

    print(f"SCORE\t{score}")
    print(f"ALNQ\t{aq}")
    print(f"ALNT\t{at}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
