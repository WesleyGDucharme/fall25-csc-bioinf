# Copyright 2008 by Bartek Wilczynski.  All rights reserved.
# Revisions copyright 2019 by Victor Lin.
# Adapted from test_Mymodule.py by Jeff Chang.
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.

"""Tests for motifs module."""

import math
from math import isfinite
import tempfile
import unittest
import os, sys

'''
# This code when uncommented runs imports the current up to date biopython module
try:
    import numpy as np
except ImportError:
    from Bio import MissingExternalDependencyError

    raise MissingExternalDependencyError(
        "Install numpy if you want to use Bio.motifs."
    ) from None

from Bio import motifs
from Bio.Seq import Seq
from Bio.motifs import thresholds as thresholds_mod
'''

# Here are the imports for the versions of the biopython motifs parts that we have locally in the repo and will be ported to codon

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
if REPO_ROOT not in sys.path:
    sys.path.append(REPO_ROOT)

from week2 import code as motifs                             # <- THIS loads week2/code/__init__.py
from week2.code import matrix as motifs_matrix
from week2.code import minimal as motifs_minimal
from week2.code import thresholds as thresholds_mod

import numpy as np
from Bio.Seq import Seq

class TestBasic(unittest.TestCase):
    """Basic motif tests."""

    def test_format(self):
        m = motifs.create([Seq("ATATA")])
        m.name = "Foo"
        s1 = format(m, "pfm")
        expected_pfm = """  1.00   0.00   1.00   0.00  1.00
  0.00   0.00   0.00   0.00  0.00
  0.00   0.00   0.00   0.00  0.00
  0.00   1.00   0.00   1.00  0.00
"""
        s2 = format(m, "jaspar")
        expected_jaspar = """>None Foo
A [  1.00   0.00   1.00   0.00   1.00]
C [  0.00   0.00   0.00   0.00   0.00]
G [  0.00   0.00   0.00   0.00   0.00]
T [  0.00   1.00   0.00   1.00   0.00]
"""
        self.assertEqual(s2, expected_jaspar)
        s3 = format(m, "transfac")
        expected_transfac = """P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        self.assertEqual(s3, expected_transfac)
        self.assertRaises(ValueError, format, m, "foo_bar")

    def test_relative_entropy(self):
        m = motifs.create([Seq("ATATA"), Seq("ATCTA"), Seq("TTGTA")])
        self.assertEqual(len(m.alignment), 3)
        self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        self.assertEqual(m.pseudocounts, {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0})
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array([1.0817041659455104, 2.0, 0.4150374992788437, 2.0, 2.0]),
            )
        )
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.8186697601117167,
                        1.7369655941662063,
                        0.5419780939258206,
                        1.7369655941662063,
                        1.7369655941662063,
                    ]
                ),
            )
        )
        m.background = None
        self.assertEqual(m.background, {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25})
        pseudocounts = math.sqrt(len(m.alignment))
        m.pseudocounts = {
            letter: m.background[letter] * pseudocounts for letter in "ACGT"
        }
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.3532586861097656,
                        0.7170228827697498,
                        0.11859369972847714,
                        0.7170228827697498,
                        0.7170228827697499,
                    ]
                ),
            )
        )
        m.background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        self.assertTrue(
            np.allclose(
                m.relative_entropy,
                np.array(
                    [
                        0.19727984803857979,
                        0.561044044698564,
                        0.20984910512125132,
                        0.561044044698564,
                        0.5610440446985638,
                    ]
                ),
            )
        )

    def test_reverse_complement(self):
        """Test if motifs can be reverse-complemented."""
        background = {"A": 0.3, "C": 0.2, "G": 0.2, "T": 0.3}
        pseudocounts = 0.5
        m = motifs.create([Seq("ATATA")])
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = format(m, "transfac")
        expected_forward = """\
P0      A      C      G      T
01      1      0      0      0      A
02      0      0      0      1      T
03      1      0      0      0      A
04      0      0      0      1      T
05      1      0      0      0      A
XX
//
"""
        self.assertEqual(received_forward, expected_forward)
        expected_forward_pwm = """\
        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.17   0.50   0.17   0.50   0.17
"""
        self.assertEqual(str(m.pwm), expected_forward_pwm)
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        expected_reverse = """\
P0      A      C      G      T
01      0      0      0      1      T
02      1      0      0      0      A
03      0      0      0      1      T
04      1      0      0      0      A
05      0      0      0      1      T
XX
//
"""
        self.assertEqual(received_reverse, expected_reverse)
        expected_reverse_pwm = """\
        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
T:   0.50   0.17   0.50   0.17   0.50
"""
        self.assertEqual(str(m.pwm), expected_reverse_pwm)
        # Same but for RNA motif.
        background_rna = {"A": 0.3, "C": 0.2, "G": 0.2, "U": 0.3}
        pseudocounts = 0.5
        m_rna = motifs.create([Seq("AUAUA")], alphabet="ACGU")
        m_rna.background = background_rna
        m_rna.pseudocounts = pseudocounts
        expected_forward_rna_counts = """\
        0      1      2      3      4
A:   1.00   0.00   1.00   0.00   1.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   0.00   1.00   0.00   1.00   0.00
"""
        self.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
        expected_forward_rna_pwm = """\
        0      1      2      3      4
A:   0.50   0.17   0.50   0.17   0.50
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.17   0.50   0.17   0.50   0.17
"""
        self.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
        expected_reverse_rna_counts = """\
        0      1      2      3      4
A:   0.00   1.00   0.00   1.00   0.00
C:   0.00   0.00   0.00   0.00   0.00
G:   0.00   0.00   0.00   0.00   0.00
U:   1.00   0.00   1.00   0.00   1.00
"""
        self.assertEqual(
            str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
        )
        expected_reverse_rna_pwm = """\
        0      1      2      3      4
A:   0.17   0.50   0.17   0.50   0.17
C:   0.17   0.17   0.17   0.17   0.17
G:   0.17   0.17   0.17   0.17   0.17
U:   0.50   0.17   0.50   0.17   0.50
"""
        self.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)
        # Same thing, but now start with a motif calculated from a count matrix
        m = motifs.create([Seq("ATATA")])
        counts = m.counts
        m = motifs.Motif(counts=counts)
        m.background = background
        m.pseudocounts = pseudocounts
        received_forward = format(m, "transfac")
        self.assertEqual(received_forward, expected_forward)
        self.assertEqual(str(m.pwm), expected_forward_pwm)
        m = m.reverse_complement()
        received_reverse = format(m, "transfac")
        self.assertEqual(received_reverse, expected_reverse)
        self.assertEqual(str(m.pwm), expected_reverse_pwm)
        # Same, but for RNA count matrix
        m_rna = motifs.create([Seq("AUAUA")], alphabet="ACGU")
        counts = m_rna.counts
        m_rna = motifs.Motif(counts=counts, alphabet="ACGU")
        m_rna.background = background_rna
        m_rna.pseudocounts = pseudocounts
        self.assertEqual(str(m_rna.counts), expected_forward_rna_counts)
        self.assertEqual(str(m_rna.pwm), expected_forward_rna_pwm)
        self.assertEqual(
            str(m_rna.reverse_complement().counts), expected_reverse_rna_counts
        )
        self.assertEqual(str(m_rna.reverse_complement().pwm), expected_reverse_rna_pwm)


    def test_minimal_meme_parser(self):
        """Parse motifs/minimal_test.meme file."""
        with open("../data/minimal_test.meme") as stream:
            record = motifs.parse(stream, "minimal")
        self.assertEqual(record.version, "4")
        self.assertEqual(record.alphabet, "ACGT")
        self.assertEqual(len(record.sequences), 0)
        self.assertEqual(record.command, "")
        self.assertEqual(len(record), 3)
        motif = record[0]
        self.assertEqual(motif.name, "KRP")
        self.assertEqual(record["KRP"], motif)
        self.assertEqual(motif.num_occurrences, 17)
        self.assertEqual(motif.length, 19)
        self.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
        self.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
        self.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
        self.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
        self.assertAlmostEqual(motif.evalue, 4.1e-09, places=10)
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertIsNone(motif.alignment)
        self.assertEqual(motif.consensus, "TGTGATCGAGGTCACACTT")
        self.assertEqual(motif.degenerate_consensus, "TGTGANNNWGNTCACAYWW")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.1684297174927525,
                        0.9432809925744818,
                        1.4307101633876265,
                        1.1549413780465179,
                        0.9308256303218774,
                        0.009164393966550805,
                        0.20124190687894253,
                        0.17618542656995528,
                        0.36777933103380855,
                        0.6635834532368525,
                        0.07729943368061855,
                        0.9838293592717438,
                        1.72489868427398,
                        0.8397561713453014,
                        1.72489868427398,
                        0.8455332015343343,
                        0.3106481207768122,
                        0.7382733641762232,
                        0.537435993300495,
                    ]
                ),
            )
        )
        self.assertEqual(motif[2:9].consensus, "TGATCGA")
        motif = record[1]
        self.assertEqual(motif.name, "IFXA")
        self.assertEqual(record["IFXA"], motif)
        self.assertEqual(motif.num_occurrences, 14)
        self.assertEqual(motif.length, 18)
        self.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
        self.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
        self.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
        self.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
        self.assertAlmostEqual(motif.evalue, 3.2e-35, places=36)
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertIsNone(motif.alignment)
        self.assertEqual(motif.consensus, "TACTGTATATATATCCAG")
        self.assertEqual(motif.degenerate_consensus, "TACTGTATATAHAWMCAG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.9632889858595118,
                        1.02677956765017,
                        2.451526420551951,
                        1.7098384161433415,
                        2.2598671267551107,
                        1.7098384161433415,
                        1.02677956765017,
                        1.391583804103081,
                        1.02677956765017,
                        1.1201961888781142,
                        0.27822438781180836,
                        0.36915366971717867,
                        1.7240522753630425,
                        0.3802185945622609,
                        0.790937683007783,
                        2.451526420551951,
                        1.7240522753630425,
                        1.3924085743645374,
                    ]
                ),
            )
        )
        self.assertEqual(motif[2:9].consensus, "CTGTATA")
        with open("../data/minimal_test.meme") as stream:
            record = motifs.parse(stream, "minimal")
        motif = record[2]
        self.assertEqual(motif.name, "IFXA_no_nsites_no_evalue")
        self.assertEqual(record["IFXA_no_nsites_no_evalue"], motif)
        self.assertEqual(motif.num_occurrences, 20)
        self.assertEqual(motif.length, 18)
        self.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
        self.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
        self.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
        self.assertAlmostEqual(motif.background["T"], 0.30569430569430567)
        self.assertAlmostEqual(motif.evalue, 0.0, places=36)
        self.assertEqual(motif.alphabet, "ACGT")
        self.assertIsNone(motif.alignment)
        self.assertEqual(motif.consensus, "TACTGTATATATATCCAG")
        self.assertEqual(motif.degenerate_consensus, "TACTGTATATAHAWMCAG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.99075309,
                        1.16078104,
                        2.45152642,
                        1.70983842,
                        2.25986713,
                        1.70983842,
                        1.16078104,
                        1.46052586,
                        1.16078104,
                        1.10213019,
                        0.29911041,
                        0.36915367,
                        1.72405228,
                        0.37696488,
                        0.85258086,
                        2.45152642,
                        1.72405228,
                        1.42793329,
                    ]
                ),
            )
        )
        self.assertEqual(motif[2:9].consensus, "CTGTATA")

    def test_meme_parser_rna(self):
        """Test if Bio.motifs can parse MEME output files using RNA."""
        with open("../data/minimal_test_rna.meme") as stream:
            record = motifs.parse(stream, "minimal")
        self.assertEqual(record.version, "4")
        self.assertEqual(record.alphabet, "ACGU")
        self.assertEqual(len(record.sequences), 0)
        self.assertEqual(record.command, "")
        self.assertEqual(len(record), 3)
        motif = record[0]
        self.assertEqual(motif.name, "KRP_fake_RNA")
        self.assertEqual(record["KRP_fake_RNA"], motif)
        self.assertEqual(motif.num_occurrences, 17)
        self.assertEqual(motif.length, 19)
        self.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
        self.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
        self.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
        self.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
        self.assertAlmostEqual(motif.evalue, 4.1e-09, places=10)
        self.assertEqual(motif.alphabet, "ACGU")
        self.assertIsNone(motif.alignment)
        self.assertEqual(motif.consensus, "UGUGAUCGAGGUCACACUU")
        self.assertEqual(motif.degenerate_consensus, "UGUGANNNWGNUCACAYWW")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        1.1684297174927525,
                        0.9432809925744818,
                        1.4307101633876265,
                        1.1549413780465179,
                        0.9308256303218774,
                        0.009164393966550805,
                        0.20124190687894253,
                        0.17618542656995528,
                        0.36777933103380855,
                        0.6635834532368525,
                        0.07729943368061855,
                        0.9838293592717438,
                        1.72489868427398,
                        0.8397561713453014,
                        1.72489868427398,
                        0.8455332015343343,
                        0.3106481207768122,
                        0.7382733641762232,
                        0.537435993300495,
                    ]
                ),
            )
        )
        self.assertEqual(motif[2:9].consensus, "UGAUCGA")
        motif = record[1]
        self.assertEqual(motif.name, "IFXA_fake_RNA")
        self.assertEqual(record["IFXA_fake_RNA"], motif)
        self.assertEqual(motif.num_occurrences, 14)
        self.assertEqual(motif.length, 18)
        self.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
        self.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
        self.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
        self.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
        self.assertAlmostEqual(motif.evalue, 3.2e-35, places=36)
        self.assertEqual(motif.alphabet, "ACGU")
        self.assertIsNone(motif.alignment)
        self.assertEqual(motif.consensus, "UACUGUAUAUAUAUCCAG")
        self.assertEqual(motif.degenerate_consensus, "UACUGUAUAUAHAWMCAG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.9632889858595118,
                        1.02677956765017,
                        2.451526420551951,
                        1.7098384161433415,
                        2.2598671267551107,
                        1.7098384161433415,
                        1.02677956765017,
                        1.391583804103081,
                        1.02677956765017,
                        1.1201961888781142,
                        0.27822438781180836,
                        0.36915366971717867,
                        1.7240522753630425,
                        0.3802185945622609,
                        0.790937683007783,
                        2.451526420551951,
                        1.7240522753630425,
                        1.3924085743645374,
                    ]
                ),
            )
        )
        self.assertEqual(motif[2:9].consensus, "CUGUAUA")

        motif = record[2]
        self.assertEqual(motif.name, "IFXA_no_nsites_no_evalue_fake_RNA")
        self.assertEqual(record["IFXA_no_nsites_no_evalue_fake_RNA"], motif)
        self.assertEqual(motif.num_occurrences, 20)
        self.assertEqual(motif.length, 18)
        self.assertAlmostEqual(motif.background["A"], 0.30269730269730266)
        self.assertAlmostEqual(motif.background["C"], 0.1828171828171828)
        self.assertAlmostEqual(motif.background["G"], 0.20879120879120877)
        self.assertAlmostEqual(motif.background["U"], 0.30569430569430567)
        self.assertAlmostEqual(motif.evalue, 0.0, places=36)
        self.assertEqual(motif.alphabet, "ACGU")
        self.assertIsNone(motif.alignment)
        self.assertEqual(motif.consensus, "UACUGUAUAUAUAUCCAG")
        self.assertEqual(motif.degenerate_consensus, "UACUGUAUAUAHAWMCAG")
        self.assertTrue(
            np.allclose(
                motif.relative_entropy,
                np.array(
                    [
                        0.99075309,
                        1.16078104,
                        2.45152642,
                        1.70983842,
                        2.25986713,
                        1.70983842,
                        1.16078104,
                        1.46052586,
                        1.16078104,
                        1.10213019,
                        0.29911041,
                        0.36915367,
                        1.72405228,
                        0.37696488,
                        0.85258086,
                        2.45152642,
                        1.72405228,
                        1.42793329,
                    ]
                ),
            )
        )
        self.assertEqual(motif[2:9].consensus, "CUGUAUA")


class MotifTestPWM(unittest.TestCase):
    """PWM motif tests."""

    with open("../data/SRF.pfm") as stream:
        m = motifs.read(stream, "pfm")

    s = Seq("ACGTGTGCGTAGTGCGT")

    def test_getitem(self):
        counts = self.m.counts
        python_integers = range(13)
        numpy_integers = np.array(python_integers)
        integers = {"python": python_integers, "numpy": numpy_integers}
        for int_type in ("python", "numpy"):
            i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12 = integers[int_type]
            msg = f"using {int_type} integers as indices"
            # slice, slice
            d = counts[i1::i2, i2:i12:i3]
            self.assertIsInstance(d, dict, msg=msg)
            self.assertEqual(len(d), 2, msg=msg)
            self.assertEqual(len(d["C"]), 4, msg=msg)
            self.assertEqual(len(d["T"]), 4, msg=msg)
            self.assertAlmostEqual(d["C"][i0], 45.0, msg=msg)
            self.assertAlmostEqual(d["C"][i1], 1.0, msg=msg)
            self.assertAlmostEqual(d["C"][i2], 0.0, msg=msg)
            self.assertAlmostEqual(d["C"][i3], 1.0, msg=msg)
            self.assertAlmostEqual(d["T"][i0], 0.0, msg=msg)
            self.assertAlmostEqual(d["T"][i1], 42.0, msg=msg)
            self.assertAlmostEqual(d["T"][i2], 3.0, msg=msg)
            self.assertAlmostEqual(d["T"][i3], 0.0, msg=msg)
            # slice, int
            d = counts[i1::i2, i4]
            self.assertIsInstance(d, dict, msg=msg)
            self.assertEqual(len(d), 2, msg=msg)
            self.assertAlmostEqual(d["C"], 1.0, msg=msg)
            self.assertAlmostEqual(d["T"], 13.0, msg=msg)
            # int, slice
            t = counts[i2, i3:i12:i2]
            self.assertIsInstance(t, tuple, msg=msg)
            self.assertAlmostEqual(t[i0], 0.0, msg=msg)
            self.assertAlmostEqual(t[i1], 0.0, msg=msg)
            self.assertAlmostEqual(t[i2], 0.0, msg=msg)
            self.assertAlmostEqual(t[i3], 0.0, msg=msg)
            self.assertAlmostEqual(t[i4], 43.0, msg=msg)
            # int, int
            v = counts[i1, i5]
            self.assertAlmostEqual(v, 1.0, msg=msg)
            # tuple, slice
            d = counts[(i0, i3), i3:i12:i2]
            self.assertIsInstance(d, dict, msg=msg)
            self.assertEqual(len(d), 2, msg=msg)
            self.assertEqual(len(d["A"]), 5, msg=msg)
            self.assertEqual(len(d["T"]), 5, msg=msg)
            self.assertAlmostEqual(d["A"][i0], 1.0, msg=msg)
            self.assertAlmostEqual(d["A"][i1], 3.0, msg=msg)
            self.assertAlmostEqual(d["A"][i2], 1.0, msg=msg)
            self.assertAlmostEqual(d["A"][i3], 15.0, msg=msg)
            self.assertAlmostEqual(d["A"][i4], 2.0, msg=msg)
            self.assertAlmostEqual(d["T"][i0], 0.0, msg=msg)
            self.assertAlmostEqual(d["T"][i1], 42.0, msg=msg)
            self.assertAlmostEqual(d["T"][i2], 45.0, msg=msg)
            self.assertAlmostEqual(d["T"][i3], 30.0, msg=msg)
            self.assertAlmostEqual(d["T"][i4], 0.0, msg=msg)
            # tuple, int
            d = counts[(i0, i3), i5]
            self.assertIsInstance(d, dict, msg=msg)
            self.assertEqual(len(d), 2, msg=msg)
            self.assertAlmostEqual(d["A"], 3.0, msg=msg)
            self.assertAlmostEqual(d["T"], 42.0, msg=msg)
            # str, slice
            t = counts["C", i2:i12:i4]
            self.assertIsInstance(t, tuple, msg=msg)
            self.assertAlmostEqual(t[i0], 45.0, msg=msg)
            self.assertAlmostEqual(t[i1], 0.0, msg=msg)
            self.assertAlmostEqual(t[i2], 0.0, msg=msg)
            # str, int
            self.assertAlmostEqual(counts["T", i4], 13.0, msg=msg)

    def test_simple(self):
        """Test if Bio.motifs PWM scoring works."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(self.s)
        self.assertEqual(6, len(result))
        # The fast C-code in Bio/motifs/_pwm.c stores all results as 32-bit
        # floats; the slower Python code in Bio/motifs/__init__.py uses 64-bit
        # doubles. The C-code and Python code results will therefore not be
        # exactly equal. Test the first 5 decimal places only to avoid either
        # the C-code or the Python code to inadvertently fail this test.
        self.assertAlmostEqual(result[0], -29.18363571, places=5)
        self.assertAlmostEqual(result[1], -38.3365097, places=5)
        self.assertAlmostEqual(result[2], -29.17756271, places=5)
        self.assertAlmostEqual(result[3], -38.04542542, places=5)
        self.assertAlmostEqual(result[4], -20.3014183, places=5)
        self.assertAlmostEqual(result[5], -25.18009186, places=5)

    def test_with_mixed_case(self):
        """Test if Bio.motifs PWM scoring works with mixed case."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(Seq("AcGTgTGCGtaGTGCGT"))
        self.assertEqual(6, len(result))
        self.assertAlmostEqual(result[0], -29.18363571, places=5)
        self.assertAlmostEqual(result[1], -38.3365097, places=5)
        self.assertAlmostEqual(result[2], -29.17756271, places=5)
        self.assertAlmostEqual(result[3], -38.04542542, places=5)
        self.assertAlmostEqual(result[4], -20.3014183, places=5)
        self.assertAlmostEqual(result[5], -25.18009186, places=5)

    def test_with_bad_char(self):
        """Test if Bio.motifs PWM scoring works with unexpected letters like N."""
        counts = self.m.counts
        pwm = counts.normalize(pseudocounts=0.25)
        pssm = pwm.log_odds()
        result = pssm.calculate(Seq("ACGTGTGCGTAGTGCGTN"))
        self.assertEqual(7, len(result))
        self.assertAlmostEqual(result[0], -29.18363571, places=5)
        self.assertAlmostEqual(result[1], -38.3365097, places=5)
        self.assertAlmostEqual(result[2], -29.17756271, places=5)
        self.assertAlmostEqual(result[3], -38.04542542, places=5)
        self.assertAlmostEqual(result[4], -20.3014183, places=5)
        self.assertAlmostEqual(result[5], -25.18009186, places=5)
        self.assertTrue(math.isnan(result[6]), f"Expected nan, not {result[6]!r}")

def _uniform_bg():
    # uniform DNA background
    return {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}

def _toy_instances():
    # Small deterministic motif, 4 columns, 4 instances:
    # pos0: A=4 ; pos1: C=4 ; pos2: G=3,C=1 ; pos3: T=3,C=1
    return [Seq("ACGT"), Seq("ACGT"), Seq("AGGT"), Seq("ACCT")]

class TestThresholds(unittest.TestCase):
    def _pssm(self):
        m = motifs.create(_toy_instances())
        # getting counts
        counts = m.counts  # FrequencyPositionMatrix
        # normalizing to PWM with explicit pseudocounts
        pwm = counts.normalize(pseudocounts=1.0)  # PositionWeightMatrix
        # converting to PSSM (log-odds) with explicit background
        return pwm.log_odds(_uniform_bg())  # PositionSpecificScoringMatrix

    def _pssm_minmax(self, pssm):
    # Portable across Biopython versions: min/max may be attributes or methods
    # Trying attributes first, then methods.
        try:
            smin = pssm.min
        except Exception:
            smin = pssm.min()
        try:
            smax = pssm.max
        except Exception:
            smax = pssm.max()
        return smin, smax

    def test_fpr_monotonic(self):
        pssm = self._pssm()
        dist = thresholds_mod.ScoreDistribution(pssm=pssm, background=_uniform_bg(), precision=1000)
        t_1 = dist.threshold_fpr(0.01)
        t_5 = dist.threshold_fpr(0.05)
        # Allowing higher FPR should LOWER the threshold
        self.assertGreater(t_1, t_5)

    def test_fnr_monotonic(self):
        pssm = self._pssm()
        dist = thresholds_mod.ScoreDistribution(pssm=pssm, background=_uniform_bg(), precision=1000)
        t_10 = dist.threshold_fnr(0.10)
        t_20 = dist.threshold_fnr(0.20)
        # Allowing higher FNR should RAISE the threshold
        self.assertLess(t_10, t_20)

    def test_balanced_threshold_brackets_and_rates(self):
        pssm = self._pssm()
        dist = thresholds_mod.ScoreDistribution(pssm=pssm, background=_uniform_bg(), precision=1000)
        t_bal = dist.threshold_balanced()
        smin, smax = self._pssm_minmax(pssm)
        # Balanced threshold should be within feasible score range
        self.assertGreaterEqual(t_bal, smin)
        self.assertLessEqual(t_bal, smax)

    def test_patser_defined_and_reasonable(self):
        pssm = self._pssm()
        dist = thresholds_mod.ScoreDistribution(pssm=pssm, background=_uniform_bg(), precision=1000)
        t_pat = dist.threshold_patser()
        self.assertTrue(isfinite(t_pat))
        smin, smax = self._pssm_minmax(pssm)
        self.assertGreaterEqual(t_pat, smin)
        self.assertLessEqual(t_pat, smax)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
