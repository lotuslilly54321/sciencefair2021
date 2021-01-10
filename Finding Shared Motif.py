from SeqIO import read_fasta_file_modif


def reader(string):
    data = read_fasta_file_modif(string)
    return data


def all_motifs(k, string):
    motifs = []
    for i in range(len(string)-k+1):
        motifs.append(string[i:i+k])
    motifs.sort()
    return motifs


def find_shared_motifs(dna_strings):
    shared = ' '
    for i in range(1, len(dna_strings[0])):
        motifs = all_motifs(len(dna_strings[0])-i, dna_strings[0])
        for motif in motifs:
            count = 0
            for dna in dna_strings:
                if motif in dna:
                    count += 1
            if count == len(dna_strings):
                if len(motif) > len(shared):
                    shared = motif
    return shared



dnas = '''
'''

data = reader('data_file')
print(data)
print('cat')
print(find_shared_motifs(data))