# import os not needed for codon


def read_fasta(path, name):
    data = []
    # building the file path manually instead of os.path.join
    file_path = path.rstrip("/\\") + "/" + name

    with open(file_path, 'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line[0] != '>':
                data.append(line)
    print(name, len(data), len(data[0]))
    # print('Sample:', data[0])
    return data


def read_data(path):
    short1 = read_fasta(path, "short_1.fasta")
    short2 = read_fasta(path, "short_2.fasta")
    long1 = read_fasta(path, "long.fasta")
    return short1, short2, long1
