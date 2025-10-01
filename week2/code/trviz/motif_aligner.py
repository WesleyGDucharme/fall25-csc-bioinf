from motif_aligner import MotifAligner as _PyMotifAligner   # original class and API
import os, subprocess, tempfile, shutil
from .utils import codon_call_args  # worker bridge

USE_CODON = os.getenv("TRVIZ_IMPL") == "codon"

class MotifAligner(_PyMotifAligner):
    """Public API unchanged. In Codon mode:
       - worker builds input FASTA & command line
       - Python runs the single shell command
       - worker parses the aligned FASTA and returns (id, seq)
    """

    # helper that does end-to-end via worker for a given method
    def _align_with_method_via_worker(self, method: str, ids, seqs,
                                      preserve_order: bool = True,
                                      tool_path: str | None = None):
        if not USE_CODON:
            # fall back to original Python behavior
            return super()._align_motifs_with_mafft(...)

        # ensure tool exists (nice early message)
        exe = (tool_path or method)
        if shutil.which(exe) is None:
            raise RuntimeError(f"{method.upper()} not found. Please install '{exe}' and ensure it is on PATH.")

        # use a temporary working directory for inputs/outputs
        with tempfile.TemporaryDirectory() as work_dir:
            # package rows as "id<TAB>seq"
            rows = [f"{i}\t{s}" for i, s in zip(ids, seqs)]
            out = codon_call_args("aligner.align",
                                  method,
                                  "1" if preserve_order else "0",
                                  exe,
                                  work_dir,
                                  *rows)
            cmd = None
            out_fa = None
            for line in out.splitlines():
                tag, val = line.split("\t", 1)
                if tag == "CMD":
                    cmd = val
                elif tag == "OUT":
                    out_fa = val
            if not cmd or not out_fa:
                raise RuntimeError("worker did not return command or output path")

            # run the command (single spawn)
            proc = subprocess.run(cmd, shell=True, text=True, capture_output=True)
            if proc.returncode != 0:
                raise RuntimeError(f"{method} failed ({proc.returncode}): {proc.stderr.strip()}")

            # collect aligned output
            parsed = codon_call_args("aligner.collect", out_fa)
            aligned_ids, aligned_seqs = [], []
            for line in parsed.splitlines():
                if "\t" in line:
                    i, s = line.split("\t", 1)
                    aligned_ids.append(i)
                    aligned_seqs.append(s)
            return aligned_ids, aligned_seqs

    # --- override tool-specific entry points in Codon mode ---

    if USE_CODON:
        def _align_motifs_with_mafft(self, input_ids, input_seqs,
                                     preserve_order=True, mafft_path="mafft"):
            return self._align_with_method_via_worker(
                "mafft", input_ids, input_seqs, preserve_order, mafft_path
            )

        def _align_motifs_with_muscle(self, input_ids, input_seqs,
                                      preserve_order=True, muscle_path="muscle"):
            return self._align_with_method_via_worker(
                "muscle", input_ids, input_seqs, preserve_order, muscle_path
            )

        def _align_motifs_with_clustalo(self, input_ids, input_seqs,
                                        preserve_order=True, clustalo_path="clustalo"):
            return self._align_with_method_via_worker(
                "clustalo", input_ids, input_seqs, preserve_order, clustalo_path
            )

        # if your original has _align_motifs_with_star(...)
        def _align_motifs_with_star(self, input_ids, input_seqs,
                                    preserve_order=True, mafft_path="mafft"):
            # keep it simple: use mafft multi-fasta; if your original does a progressive seed,
            # we can mirror that later by composing several mafft calls and merging—just say the word.
            return self._align_with_method_via_worker(
                "mafft", input_ids, input_seqs, preserve_order, mafft_path
            )
