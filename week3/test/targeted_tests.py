from pathlib import Path
import numpy as np
import pytest

# choose implementation via env var: PHYLO_IMPL=biotite|local
import os
if os.getenv("PHYLO_IMPL", "biotite").lower() == "local":
    import phylo_local_files as phylo
else:
    import biotite.sequence.phylo as phylo

# data helpers: look in week3/data
BASE_DIR = Path(__file__).resolve().parents[1]
SEQ_DIR = BASE_DIR / "data"

@pytest.fixture
def distances():
    """
    Distances are based on "Dendrogram of the BLOSUM62 matrix"
    with the small modification M[i,j] += i+j to reduce ambiguity.
    """
    return np.loadtxt(SEQ_DIR / "distances.txt", dtype=int)

@pytest.fixture
def upgma_newick():
    with open(SEQ_DIR / "newick_upgma.txt", "r", encoding="utf-8") as f:
        return f.read().strip()

@pytest.fixture
def tree(distances):
    return phylo.upgma(distances)

def test_upgma(tree, upgma_newick):
    """Compare the results of `upgma()` with DendroUPGMA."""
    ref_tree = phylo.Tree.from_newick(upgma_newick)
    # Compare pairwise distances & topology with tolerance for FP rounding
    for i in range(len(tree)):
        for j in range(len(tree)):
            assert tree.get_distance(i, j) == pytest.approx(
                ref_tree.get_distance(i, j), abs=1e-3
            )
            assert tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                i, j, topological=True
            )

def test_neighbor_joining():
    """Compare `neighbor_joining()` with a known tree."""
    dist = np.array([
        [0, 5, 4, 7, 6, 8],
        [5, 0, 7, 10, 9, 11],
        [4, 7, 0, 7, 6, 8],
        [7, 10, 7, 0, 5, 9],
        [6, 9, 6, 5, 0, 8],
        [8, 11, 8, 9, 8, 0],
    ])

    ref_tree = phylo.Tree(
        phylo.TreeNode(
            [
                phylo.TreeNode(
                    [
                        phylo.TreeNode(
                            [phylo.TreeNode(index=0), phylo.TreeNode(index=1)],
                            [1, 4],
                        ),
                        phylo.TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                phylo.TreeNode(
                    [phylo.TreeNode(index=3), phylo.TreeNode(index=4)],
                    [3, 2],
                ),
                phylo.TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )

    test_tree = phylo.neighbor_joining(dist)
    assert test_tree == ref_tree

def test_distances(tree):
    """Sanity checks on distances in the UPGMA tree."""
    # Equal root distance for all leaves (UPGMA property)
    dist0 = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist0
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10
