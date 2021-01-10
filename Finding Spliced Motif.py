from SeqIO import read_fasta_file_modif


def reader(string):
    data = read_fasta_file_modif(string)
    return data


def spliced_motif(s, t):
    found = ''
    to_find = t
    indexes = []
    for i in range(len(s)):
        if len(to_find) > 0:
            if s[i] == to_find[0]:
                indexes.append(i+1)
                to_find = to_find[1:]
                print to_find
    if len(indexes) != len(t):
        return False
    return True


data = None#reader('data_file')
dna = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAA'#data[0]
motif = "CCCC"#data[1]
print dna
print motif
print spliced_motif(dna, motif)