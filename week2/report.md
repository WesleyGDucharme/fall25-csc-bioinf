Overall the TRViz codebase was ported to run with a Codon backend while keeping the original Python API and tests intact. The final setup is more of a hybrid that uses the code I ported to Codon in tandem with the python code:  

- The original Python/Cython behavior remains the reference.  
- When TRVIZ_IMPL=codon, Python routes certain operations to a persistent Codon worker (trviz_codon_worker) for speed and lower overhead.  
- A Codon-only CLI (trviz_codon_cli) mirrors trviz/main.py.  
  
All official unit tests (test_decomposer.py, test_refinement.py) pass in both Python and Codon modes. We also added small Codon smoke tests for the worker/CLI that specifically targets Codon ported code and runs it separately.  
One important caveat: I included a Codon “DP” fallback for decomposition that is greedy (by design) and does not perfectly match the Cython DP. In our Python path we still prefer Cython (when available) to preserve exact behavior.  

Steps  
-----
1) Cloning & Baseline  

Cloned the week-2 repo, set up a venv, verified the provided Python tests pass locally.  

Confirmed Cython-backed DP path runs as expected when available; else it falls back to Python DP.  

2) Design & Bridge  

Decided on a persistent worker model for Codon (single process, line protocol over stdin/stdout).  

Wrote trviz/_codon_bridge.py and exposed codon_call_args(cmd, *args) to Python callers.  

This avoided per-call Codon startup overhead and enabled “swappable” backends via TRVIZ_IMPL.  

3) Porting Codon  

trviz_codon/utils.codon:  
FASTA parsing (no BioPython), Levenshtein, distance/score matrices, add_padding, sorting, motif marks, etc.  

trviz_codon/motif_encoder.codon:  
Private motif thresholding & encoding; mirrors Python output format.  

trviz_codon/motif_aligner.codon:  
Minimal MAFFT wrapper (writes FASTA, runs /usr/bin/mafft, reads output).  

trviz_codon/decomposer.codon:  
Greedy DP fallback stub for decomposer.dp to keep the port complete.  

trviz_codon/worker.codon:  
Exposes commands (utils.get_fasta, utils.levenshtein, encoder.encode, decomposer.refine, decomposer.dp, etc.) using a small, consistent protocol.  

trviz_codon/main.codon:  
Full CLI equivalent of trviz/main.py (arg parsing, FASTA IO, optional MAFFT, write TSV).  

4) Python Shims & Routing  

trviz/utils.py: uses Codon worker for FASTA & helpers when TRVIZ_IMPL=codon, else pure Python.  

trviz/motif_encoder.py: calls Codon worker for threshold/encode; computes score matrix same as Python.  

trviz/decomposer.py:  

Keeps DP_CY (Cython) as primary if requested/available.  

Adds Codon-backed refine() and Codon “DP” fallback (used only if asked; otherwise Cython/Python DP).  

trviz/motif_aligner.py: optional path to Codon MAFFT for smoke tests.  

5) Testing & Automation  

Kept original tests and made sure they pass in both modes, running from week2/test when going through the porting process:  
PYTHONPATH=../code TRVIZ_IMPL=python pytest -q test_decomposer.py test_refinement.py  
PYTHONPATH=../code TRVIZ_IMPL=codon  pytest -q test_decomposer.py test_refinement.py  

Added test/test_codon_port_smoke.py to exercise the worker/CLI specific pieces (FASTA, motif counter, score matrix, encoder/threshold, refine).  

Wrote week2/evaluate.sh to:  
    
- Build the worker & CLI  

- Run tests in Python mode  

- Run tests in Codon mode  

- Run the smoke tests  

- Run a CLI demo and print the first lines of the TSV  

6) CI (GitHub Actions)  

Workflow installs Codon into $HOME/.codon, sets CODON_PYTHON, and runs week2/evaluate.sh on push.  
Ensured paths use $HOME (shell env) and ${{ env.CODON_PYTHON }} (Actions env) correctly.  

Gotchas  
-------

Codon typing rules that I ran into alot.  
  
- No calls in type annotations error popped up in early versions I wrote (e.g., use Tuple[str,str] not -> (str, str)).  

- Dict.get(key, default) (two-arg) isn’t available; used v = d.get(key); v = default if v is None else v.  

- List[Optional[str]] vs List[str] mismatches surface when parsing files; so added a coerce_pairs() helper.  

Protocol edge cases  

- Worker must always emit the OK\n<count>\n... header. A missing newline caused early failures.  

- Normalized outputs so Python smoke tests expect exact line endings ("1\n" vs "1").  

I/O differences  

- Codon File doesn’t have .readline(); iterate with for line in fh.  

- When reading/writing FASTA in Codon, this avoids BioPython entirely.  
  
Imports & stdlib  
  
- Use import os and os.system(cmd) in Codon for shelling out   

Decomposer DP parity  

- The Codon decomposer.dp is greedy and not a full DP. It’s included for portability/completeness.  
  - This ported greedy version does not give the exact output as Cython/Python fallback version does but it is very close (I think based on looking at its outputs it is often off by 1 letter) and I tried to implement it in an attempt to replicate their results using just codon.  

- For exact behavior, keep using the Cython DP (DP_CY mode) in Python. Which is what we do in the main running of the tests in evaluate.sh.  


File structure and contents outline  
-----------------------------------

week2/  
├─ code/  
│  ├─ trviz/                     # Python API (unchanged entry points)  
│  │  ├─ _codon_bridge.py        # persistent worker bridge (always importable)  
│  │  ├─ utils.py                # routes to Codon when TRVIZ_IMPL=codon  
│  │  ├─ decomposer.py           # Cython DP preferred; Codon refine + DP fallback  
│  │  ├─ motif_encoder.py        # Codon for threshold/encode  
│  │  ├─ motif_aligner.py        # optional Codon MAFFT  
│  │  └─ main.py                 # original Python CLI  
│  ├─ trviz_codon/               # Codon sources  
│  │  ├─ utils.codon  
│  │  ├─ motif_encoder.codon  
│  │  ├─ motif_aligner.codon  
│  │  ├─ decomposer.codon        # greedy DP fallback  
│  │  ├─ worker.codon            # persistent worker  
│  │  ├─ main.codon              # standalone Codon CLI  
│  │  └─ __init__.py / __init__.codon  
│  ├─ trviz_codon_worker         # built worker binary (from evaluate.sh)  
│  └─ trviz_codon_cli            # built CLI binary (from evaluate.sh)  
├─ data/  
│  └─ demo.fasta  
├─ test/  
│  ├─ test_decomposer.py  
│  ├─ test_refinement.py  
│  └─ test_codon_port_smoke.py   # added smoke tests for the Codon port components  
└─ evaluate.sh                    # builds worker & CLI; runs original tests and the new small smoke tests  

The other .py files you see in week2/code are the original .py files from Trviz.  

Estimated time to complete: 
Spent more than 20 hrs on BioPython port before I gave up on that one and swapped to Trviz, where I spent ~ 18 hours on this one.
