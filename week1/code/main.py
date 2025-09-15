from dbg import DBG
from utils import read_data
import sys
# import os do not need this import anymore

# sys.setrecursionlimit(1000000) removing this line for Codon conversion as Codon does not use CPython's recursion guard.

# Raising recursion limit only if the runtime supports it (CPython does; Codon ignores)
if hasattr(sys, "setrecursionlimit"):
    try:
        sys.setrecursionlimit(1_000_000)
    except Exception:
        pass

if __name__ == "__main__":
    argv = sys.argv
    # short1, short2, long1 = read_data(argv[1])
    #with open(argv[1] + '/contig.fasta', 'w') as f:
    data_dir = argv[1].rstrip("/\\") #nomralize trailing slash
    short1, short2, long1 = read_data(data_dir)

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])
    # dbg.show_count_distribution()
    out_path = data_dir + "/contig.fasta" # avoiding os.path.join for Codon
    with open(out_path,  'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            print(i, len(c))
            f.write(f'>contig_{i}\n')
            f.write(c + '\n')
