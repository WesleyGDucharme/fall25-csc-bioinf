Week 3 Report — Porting biotite.sequence.phylo (subset) to Codon  
----------------------------------------------------------------
Ported just enough of Biotite’s phylo functionality to Codon to pass three targeted tests:  

- test_distances  
- test_upgma  
- test_neighbor_joining  

Kept a clean separation between Python baseline (Biotite) and the Codon port, and made evaluate.sh that runs both and prints the table showign the runtimes.  
 
Steps  
-----  

Baseline & isolation  
- Isolated the three required tests into test/targeted_tests.py and wired data files to week3/data/.  
- Verified the Python baseline using Biotite: python test/run_python_tests.py.  

Codon milestones (build up minimal API)  
 
Tree/TreeNode core (distances, equality)  
- Implemented TreeNode immutability constraints, parent/child wiring, and distance paths.  
- Implemented equality via a canonical, order-insensitive structural key (including edge lengths).  
- Implemented topological (edge-count) distance with a dedicated helper to avoid float/int issues.  
- Implemented Tree.from_newick() and to_newick() (subset sufficient for tests).  

Neighbor-Joining (NJ)  
- Implemented classical NJ: compute r[i], minimize Q, compute limb lengths, update distance matrix, merge clusters.  
- Hand-built reference tree matches expected topology for the 6×6 case.  

UPGMA  
- Implemented size-weighted merges and correct parent heights (h_parent = 0.5 * D[i][j]), with child edges Li = h_parent - height[i].  
- Verified pairwise distances vs. reference Newick from newick_upgma.txt.  

Codon test harness  
- Rewrote the three tests as code/test_phylo.codon, reusing the data files.  
- Added code/run_codon_tests.codon to time all three tests and print either a table or just ms (--just-ms).  

Design Choices (Often resulted in solutions to errors I was stuck on)
---------------------------------------------------------------------

- One return type per function. Ensured functions like get_distance() always return float; tests cast to int when comparing topological values.  
- Dict keys for nodes. Codon requires both __hash__ and __eq__; used identity semantics (hash by id(self), equality by self is other) for path maps.  
- Type strictness: When in doubt, casted dict lookups to float at the use site to satisfy Codon’s realization (float(a_len[cur])).  
- No closures for parser state. The Newick parser was written as a small class with explicit fields (text, pos) to avoid capture issues.  
- File paths & CWD. Codon tests rely on running inside week3/code/ (or we resolve paths carefully). The evaluate.sh runs the Codon step from that directory to avoid path ambiguity.  

Repo layout  
-----------  

week3/  
├─ code/  
│  ├─ tree.codon          # Tree/TreeNode + Newick + distances/equality  
│  ├─ upgma.codon         # UPGMA (weighted cluster merge)  
│  ├─ nj.codon            # Neighbor-Joining  
│  ├─ phylo.codon         # aggregator: exports Tree, TreeNode, upgma, neighbor_joining  
│  ├─ test_phylo.codon    # Codon reimplementation of the 3 tests  
│  ├─ run_codon_tests.codon  # times all 3 Codon tests; --just-ms prints only ms  
│  └─ (Cython sources kept for reference: .pyx files)  
├─ data/  
│  ├─ distances.txt  
│  └─ newick_upgma.txt  
├─ test/  
│  ├─ targeted_tests.py       # Python tests with env-switchable impl (biotite | local)  
│  ├─ run_python_tests.py     # times the Python tests; --just-ms prints only ms  
│  └─ phylo_local_files.py    # (experimental only) local Cython import path (not used by CI and doesn't work)  
└─ evaluate.sh                # runs both timing scripts and prints the table  

Estimated time to complete  
--------------------------  
Around 15 hours.  
