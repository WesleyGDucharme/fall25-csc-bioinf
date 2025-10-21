Week 4 Report — Sequence Alignment in Codon & Python  
----------------------------------------------------  

Implemented four alignment algorithms in Codon with identical versions in Python, plus a shared test harness and CI timing script.  

Algorithms (scoring fixed unless noted):  
- Global (Needleman–Wunsch, linear gaps) — match +3, mismatch −3, gap −2.  
- Local (Smith–Waterman, linear gaps) — same linear scores, with local reset to 0.  
- Global (Gotoh, affine gaps) — gap open −5, gap extend −1.  
- Semi-global / Fitting (affine gaps) — query must fully align, target prefixes/suffixes are free; open −5, extend −1.  

Determinism across all methods:  
- Tie-break order: diag > up > left (strict > so the first wins).  
- Local start cell: first (row-major) global max.  
- Fitting terminal cell: best in the last row, leftmost on ties.  

Parities verified for q1..q5 vs t1..t5 (multi-FASTA) and for MT_human/MT_orang in windows/truncations. CI prints a timing table for both languages.  

Steps  
-----  

1) Codon implementations  
Wrote global_align_linear, local_align_linear, global_align_affine (Gotoh with M/X/Y), and fitting_align_affine (semi-global with affine).  
Removed tuple return annotations (Codon: “cannot use calls in type signatures”); let Codon infer.  
Avoided List[str] = [] on locals; used aq = []/at = [] so Codon infers element types from append.  

2) Python parity  
Mirrored the Codon logic line-for-line (constants, tie-breaks, inits, traceback).  
Verified with diff -u on all q/t pairs and on multiple MT windows/slices.  

3) Test harnesses (shared behavior)  
Input modes: raw strings or FASTA (--fasta).  
Multi-FASTA selection: --index N (1-based qN/tN), or --qid qK + --tid tK.  
MT windows / memory control: --limit N, --slice-q a:b, --slice-t a:b (1-based inclusive).  

Output (3 lines):  
SCORE  <int>  
ALNQ   <aligned_query>  
ALNT   <aligned_target>  

4) Evaluate script & timing  
Prebuild Codon once: codon build -release -o .bin/codon_test code/codon_test.codon (avoids ~2s compile per call).  
Robust timing in CI: /usr/bin/time -f %e → ms; fallback to date if needed.  
Clear labels: globalL, local, fittingA, globalA.  
Includes asymmetric fittingA-mt_orang (orang→human) since fitting is directional.  

5) Layout & data  
week4/  
├─ code/  
│  ├─ codon_align.codon      # 4 Codon algorithms  
│  ├─ codon_test.codon       # Codon CLI harness (multi-FASTA, slices/limits)  
│  └─ py_align.py            # 4 Python algorithms (parity with Codon)  
├─ data/  
│  ├─ MT-human.fa            # full mtDNA  
│  ├─ MT-orang.fa            # full mtDNA  
│  ├─ q1.fa                  # contains q1..q5  
│  └─ t1.fa                  # contains t1..t5  
├─ test/  
│  └─ py_test.py             # Python CLI harness (same flags & output)  
└─ evaluate.sh               # builds Codon binary, times Codon & Python, prints table  
