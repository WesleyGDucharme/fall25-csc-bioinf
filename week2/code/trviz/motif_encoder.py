from motif_encoder import *


from motif_encoder import MotifEncoder as _PyMotifEncoder  # the original Python class
import os
from .utils import INDEX_TO_CHR, PRIVATE_MOTIF_LABEL, get_score_matrix
from .utils import codon_call_args

USE_CODON = os.getenv("TRVIZ_IMPL") == "codon"

class MotifEncoder(_PyMotifEncoder):
    """Thin wrapper that keeps the original Python behavior, but overrides
    selected methods when TRVIZ_IMPL=codon."""
    if USE_CODON:
        def find_private_motif_threshold(self, decomposed_vntrs, label_count=None, auto=False):
            rows = [",".join(v) for v in decomposed_vntrs]
            lc = str(label_count if (label_count is not None) else -1)
            out = codon_call_args("encoder.find_threshold", lc, *rows)
            for line in out.splitlines():
                k, v = line.split("\t", 1)
                if k == "THRESH":
                    return int(v)
            raise ValueError("threshold not returned by worker")

        def encode(self, decomposed_vntrs, label_count=None, auto=False,
                   match_score=2,
                   mismatch_score_for_edit_dist_of_1=-1,
                   mismatch_score_for_edit_dist_greater_than_1=-2,
                   gap_open_penalty=1.5,
                   gap_extension_penalty=0.6,
                   motif_map_output_file=None):
            rows = [",".join(v) for v in decomposed_vntrs]
            lc = str(label_count if (label_count is not None) else -1)
            af = "1" if auto else "0"
            out = codon_call_args("encoder.encode", lc, af, *rows)

            motif_to_symbol = {}
            symbol_to_motif = {}
            motif_counter = {}
            encoded = []
            private_thresh = None

            for line in out.splitlines():
                parts = line.split("\t")
                tag = parts[0]
                if tag == "THRESH":
                    private_thresh = int(parts[1])
                elif tag == "MAP":
                    mot, sym, cnt = parts[1], parts[2], int(parts[3])
                    motif_to_symbol[mot] = sym
                    symbol_to_motif[sym] = mot
                    motif_counter[mot] = cnt
                elif tag == "ENC":
                    encoded.append(parts[1])

            if private_thresh is None:
                raise ValueError("encode: worker did not return THRESH")

            score_matrix = get_score_matrix(
                symbol_to_motif,
                match_score=match_score,
                mismatch_score_for_edit_dist_of_1=mismatch_score_for_edit_dist_of_1,
                mismatch_score_for_edit_dist_greater_than_1=mismatch_score_for_edit_dist_greater_than_1,
                gap_open_penalty=gap_open_penalty,
                gap_extension_penalty=gap_extension_penalty,
            )

            # set the same attributes the Python version sets
            self.motif_to_symbol = motif_to_symbol
            self.symbol_to_motif = symbol_to_motif
            self.motif_counter = motif_counter
            self.score_matrix = score_matrix
            self.private_motif_threshold = private_thresh

            if motif_map_output_file:
                with open(motif_map_output_file, "w") as f:
                    for m, s in motif_to_symbol.items():
                        f.write(f"{m}\t{s}\t{motif_counter.get(m, 0)}\n")

            return encoded
