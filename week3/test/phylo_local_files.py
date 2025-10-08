from pathlib import Path
import sys, os

# --- ensure distutils import works on Python 3.12 ---
try:
    import distutils  # noqa: F401
except ModuleNotFoundError:
    import setuptools
    sys.modules['distutils'] = setuptools._distutils  # provide shim
# (optional but helps some environments)
os.environ.setdefault("SETUPTOOLS_USE_DISTUTILS", "local")

# Makes week3/code importable
SRC = Path(__file__).resolve().parents[1] / "code"
sys.path.insert(0, str(SRC))

# Compile & import .pyx on the fly
import numpy as np
import pyximport
pyximport.install(language_level=3, setup_args={"include_dirs": np.get_include()})

# Import Cython modules
import tree as _tree
import upgma as _upgma
import nj as _nj

# Re-exports a Biotite-like API so tests can do `import phylo_local as phylo`
Tree = _tree.Tree
TreeNode = _tree.TreeNode
TreeError = getattr(_tree, "TreeError", Exception)

upgma = _upgma.upgma
neighbor_joining = _nj.neighbor_joining

__all__ = ["Tree", "TreeNode", "TreeError", "upgma", "neighbor_joining"]
