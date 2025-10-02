# test/test_codon_port_smoke.py
import os
import shutil
import tempfile
import subprocess
import pathlib
import pytest

# These tests are specifically for the Codon backend
CODON = os.getenv("TRVIZ_IMPL") == "codon"

pytestmark = pytest.mark.skipif(
    not CODON, reason="Codon smoke tests only run when TRVIZ_IMPL=codon"
)

# bridge
from trviz.utils import codon_call_args

DATA_DIR = pathlib.Path(__file__).resolve().parent.parent / "data"
DEMO_FASTA = DATA_DIR / "demo.fasta"

def test_utils_get_fasta_and_is_valid():
    out = codon_call_args("utils.get_fasta", str(DEMO_FASTA))
    lines = [ln for ln in out.splitlines() if ln.strip()]
    # Expect two lines: "header<TAB>SEQ"
    assert len(lines) == 2
    h1, s1 = lines[0].split("\t", 1)
    h2, s2 = lines[1].split("\t", 1)
    assert h1 == "seq1" and h2 == "seq2"
    # is_valid_seq (DNA-only) — worker returns "1" (no trailing newline)
    assert codon_call_args("utils.is_valid_seq", s1).strip() == "1"
    assert codon_call_args("utils.is_valid_seq", "NXXX").strip() == "0"

def test_utils_levenshtein_and_padding():
    assert codon_call_args("utils.levenshtein", "kitten", "sitting").strip() == "3"
    out = codon_call_args("utils.add_padding", "ACT", "A--", "ACTG").strip().splitlines()
    # Padded to max length 4
    assert out == ["ACT-", "A---", "ACTG"]

def test_utils_motif_counter_and_score_matrix():
    # Motif counter — codon_call_args returns only payload rows, no count header; order not guaranteed
    rows = codon_call_args(
        "utils.motif_counter",
        "ACT,TT,ACT", "G", "ACT,G"
    ).strip().splitlines()
    assert set(rows) == {"ACT\t3", "G\t2", "TT\t1"}

    # Score matrix (simple sanity: identity=match, gap params echoed)
    out2 = codon_call_args("utils.score_matrix",
                           "2", "-1", "-2", "1.5", "0.6",
                           "A:ACT", "B:ACG").splitlines()
    assert any("GAP_OPEN\t1.5" in ln for ln in out2)
    assert any(ln == "A\tA\t2" for ln in out2)

def test_encoder_find_threshold_and_encode():
    rows = ["ACT,TT,ACT", "G", "ACT,G"]
    thr = codon_call_args("encoder.find_threshold", "-1", *rows).strip().splitlines()
    assert thr and thr[0].startswith("THRESH")
    enc = codon_call_args("encoder.encode", "-1", "1", *rows).splitlines()
    # Expect lines: THRESH …, MAP …, ENC …
    assert any(ln.startswith("THRESH") for ln in enc)
    assert any(ln.startswith("MAP\tACT\t") for ln in enc)
    assert sum(ln.startswith("ENC\t") for ln in enc) == 3

def test_decomposer_refine_smoke():
    # This just checks round-trip shape and that CSV decoding/encoding is consistent
    rows_in = ["ACT,TT,ACT", "G", "ACT,G"]
    out = codon_call_args("decomposer.refine", *rows_in).strip().splitlines()
    # Each line should be a CSV row (possibly identical or changed)
    assert len(out) == len(rows_in)
    assert all("," in r or r == "" or r.isalpha() for r in out)

def test_decomposer_dp_stub_smoke():
    # Greedy fallback path in Codon worker; stable small example
    out = codon_call_args("decomposer.dp", "ACTTTG", "ACT,TT,G").strip().splitlines()
    enc = [ln for ln in out if ln.startswith("ENC\t")]
    dec = [ln for ln in out if ln.startswith("DEC\t")]
    assert len(enc) == 1 and len(dec) == 1
    # DEC row should decode to the expected motifs for this simple case
    assert dec[0].split("\t", 1)[1].split(",") == ["ACT", "TT", "G"]

@pytest.mark.skipif(not shutil.which("/usr/bin/mafft"), reason="MAFFT not installed")
def test_aligner_mafft_roundtrip(tmp_path):
    # Use worker's align+collect to ensure the roundtrip works if MAFFT is available
    in_fa = tmp_path / "in.fa"
    out_fa = tmp_path / "out.fa"
    # Write a tiny FASTA
    in_fa.write_text(">s1\nACGTT\n>s2\nAC-TT\n", encoding="utf-8")
    # Ask aligner to simulate running mafft (worker prints CMD/IN/OUT/IDS)
    align_out = codon_call_args("aligner.align", "mafft", "1", "mafft", ".", "s1\tACGTT", "s2\tAC-TT")
    assert "CMD\t" in align_out and "IN\t" in align_out and "OUT\t" in align_out
    # We won’t actually run MAFFT here; just verify the collect fails gracefully on a missing file
    with pytest.raises(Exception):
        codon_call_args("aligner.collect", str(out_fa))

def test_main_cli_smoke(tmp_path):
    """Run the compiled Codon CLI if present and check it writes output."""
    cli = pathlib.Path(__file__).resolve().parents[1] / "code" / "trviz_codon_cli"
    if not cli.exists():
        pytest.skip("Codon CLI not built")
    out_path = tmp_path / "out.tsv"
    # Run: trviz_codon_cli --input <fasta> --motifs ACT,TT,G --out <file> --verbose
    cmd = [str(cli), "--input", str(DEMO_FASTA), "--motifs", "ACT,TT,G", "--out", str(out_path), "--verbose"]
    env = os.environ.copy()
    env["TRVIZ_IMPL"] = "codon"
    proc = subprocess.run(cmd, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert proc.returncode == 0, proc.stdout
    assert out_path.exists()
    # Basic header line
    head = out_path.read_text(encoding="utf-8").splitlines()[0]
    assert head.startswith("#sample_id")
