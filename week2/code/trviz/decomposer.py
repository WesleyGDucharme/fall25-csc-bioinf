# Keeping public API the same, swaping refine() to Codon when TRVIZ_IMPL=codon.

import os
from decomposer import Decomposer as _PyDecomposer  # the original class in ../decomposer.py

USE_CODON = os.getenv("TRVIZ_IMPL") == "codon"

class Decomposer(_PyDecomposer):
    # Inherit everything; override refine as a staticmethod when in Codon mode
    if USE_CODON:
        @staticmethod
        def refine(decomposed_trs, verbose: bool = False):
            # rows -> "m1,m2,..."
            from .utils import codon_call_args
            rows = [",".join(tr) for tr in decomposed_trs]
            out = codon_call_args("decomposer.refine", *rows)
            # Worker returns:
            # OK\n<cnt>\nrow\nrow\n...
            # codon_call_args already strips the OK header and count and returns just payload lines
            refined = []
            for line in out.splitlines():
                refined.append(line.split(",") if line else [])
            return refined
