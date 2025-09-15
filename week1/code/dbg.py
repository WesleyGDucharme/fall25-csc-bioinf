import copy
from typing import Dict, Set, List, Optional

# --- CPython set interop (Codon) / fallback (Python) ---
try:
    import python  # enable CPython interop in Codon

    @python
    def _py_set_new():
        return set()

    @python
    def _py_set_add(s, x):
        s.add(x)
        return None

    @python
    def _py_set_list(s):
        return list(s)

    @python
    def _py_set_sub(a, b):
        return a - b

    # >>> Typed wrappers so Codon sees Set[int] and List[int] <<<
    def py_set_new_int() -> Set[int]:
        s: Set[int] = _py_set_new()
        return s

    def py_set_add_int(s: Set[int], x: int) -> None:
        _py_set_add(s, x)

    def py_set_list_int(s: Set[int]) -> List[int]:
        lst: List[int] = _py_set_list(s)
        return lst

    def py_set_sub_int(a: Set[int], b: Set[int]) -> Set[int]:
        r: Set[int] = _py_set_sub(a, b)
        return r

except Exception:
    # running under plain Python — just use built-ins with matching types
    def py_set_new_int() -> Set[int]:
        return set()

    def py_set_add_int(s: Set[int], x: int) -> None:
        s.add(x)

    def py_set_list_int(s: Set[int]) -> List[int]:
        return list(s)

    def py_set_sub_int(a: Set[int], b: Set[int]) -> Set[int]:
        return a - b

def reverse_complement(key):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    key = list(key[::-1])
    for i in range(len(key)):
        key[i] = complement[key[i]]
    return ''.join(key)


class Node:
    _children: Set[int]
    _count: int
    kmer: str
    visited: bool
    depth: int
    max_depth_child: Optional[int]  # None sentinel like the original

    def __init__(self, kmer):
        self._children = py_set_new_int()   # CPython set under Codon
        self._count = 0
        self.kmer = kmer
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def add_child(self, kmer: int):
        py_set_add_int(self._children, kmer)

    def increase(self):
        self._count += 1

    def reset(self):
        self.visited = False
        self.depth = 0
        self.max_depth_child = None

    def get_count(self) -> int:
        return self._count

    def get_children(self) -> List[int]:
        return py_set_list_int(self._children)

    def remove_children(self, target: Set[int]):
        self._children = py_set_sub_int(self._children, target)


class DBG:
    def __init__(self, k, data_list):
        self.k = k
        self.nodes: Dict[int, Node] = {}
        # private
        self.kmer2idx: Dict[str, int] = {}
        self.kmer_count = 0
        self.order: List[int] = []  # explicit root insertion order
        # build
        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list):
        # check data list
        assert len(data_list) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list):
        for data in data_list:
            for original in data:
                rc = reverse_complement(original)
                for i in range(len(original) - self.k - 1):
                    self._add_arc(original[i: i + self.k], original[i + 1: i + 1 + self.k])
                    self._add_arc(rc[i: i + self.k], rc[i + 1: i + 1 + self.k])

    def show_count_distribution(self):
        count = [0] * 30
        for idx in self.nodes:
            count[self.nodes[idx].get_count()] += 1
        print(count[0:10])
        # plt.plot(count)
        # plt.show()

    def _add_node(self, kmer):
        # match original behavior: no k-mer length guard here
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.order.append(self.kmer_count)  # track insertion order for roots
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx

    def _add_arc(self, kmer1, kmer2):
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)

    def _get_count(self, child):
        return self.nodes[child].get_count()

    def _get_sorted_children(self, idx):
        children = self.nodes[idx].get_children()
        children.sort(key=self._get_count, reverse=True)  # stable; ties follow set iteration
        return children

    def _get_depth(self, idx):
        if not self.nodes[idx].visited:
            self.nodes[idx].visited = True
            children = self._get_sorted_children(idx)
            max_depth, max_child = 0, None
            for child in children:
                depth = self._get_depth(child)
                if depth > max_depth:
                    max_depth, max_child = depth, child
            self.nodes[idx].depth, self.nodes[idx].max_depth_child = max_depth + 1, max_child
        return self.nodes[idx].depth

    def _reset(self):
        for idx in self.nodes.keys():
            self.nodes[idx].reset()

    def _get_longest_path(self) -> List[int]:
        max_depth, max_idx = 0, None
        for idx in self.order:  # explicit insertion order for roots
            depth = self._get_depth(idx)
            if depth > max_depth:
                max_depth, max_idx = depth, idx

        path: List[int] = []
        while max_idx is not None:
            path.append(max_idx)
            max_idx = self.nodes[max_idx].max_depth_child
        return path

    def _delete_path(self, path: List[int]) -> None:
        path_set: Set[int] = set(path)
        for idx in path:
            del self.nodes[idx]
        # keep the explicit insertion order in sync
        self.order = [i for i in self.order if i not in path_set]
        for idx in list(self.nodes.keys()):
            self.nodes[idx].remove_children(path_set)

    def _concat_path(self, path):
        if len(path) < 1:
            return None
        concat = copy.copy(self.nodes[path[0]].kmer)
        for i in range(1, len(path)):
            concat += self.nodes[path[i]].kmer[-1]
        return concat

    def get_longest_contig(self):
        # reset params in nodes for getting longest path
        self._reset()
        path = self._get_longest_path()
        contig = self._concat_path(path)
        self._delete_path(path)
        return contig
