from SeqIO import read_fasta_file_modif
from itertools import combinations


def reader(string):
    data = read_fasta_file_modif(string)
    return data


def rSubset(arr, r):
    # return list of all subsets of length r
    # to deal with duplicate subsets use
    # set(list(combinations(arr, r)))
    results = set(list(combinations(arr, r)))
    strings = []
    for result in results:
        string = ''
        for r in result:
            string += r
        strings.append(string)
    strings.sort()
    return strings


def all_motifs(k, string):
    motifs = []
    for i in range(len(string)-k+1):
        motifs.append(string[i:i+k])
    motifs.sort()
    return motifs


def find_shared_motifs(dna_strings):
    shared_motifs = []
    for i in range(1, len(dna_strings[0])):
        motifs = all_motifs(len(dna_strings[0])-i, dna_strings[0])
        for motif in motifs:
            count = 0
            for dna in dna_strings:
                if motif in dna:
                    count += 1
            if count == len(dna_strings):
                shared_motifs.append(motif)
    return shared_motifs


def spliced_motif(s, t):
    to_find = t
    indexes = []
    for i in range(len(s)):
        if len(to_find) > 0:
            if s[i] == to_find[0]:
                indexes.append(i+1)
                to_find = to_find[1:]
    if len(indexes) != len(t):
        return False
    return True


def print_array(whole):
    for array in whole:
        print array


def shared_spliced_array(t1, t2):
    whole_array = []
    for i in range(len(t1)+1):
        lst = []
        for j in range(len(t2)+1):
            lst.append(0)
        whole_array.append(lst)
    for i in range(1, len(t1)+1):
        for j in range(1, len(t2)+1):
            if t1[i-1] == t2[j-1]:
                whole_array[i][j] = whole_array[i-1][j-1] + 1
            else:
                whole_array[i][j] = max(whole_array[i-1][j],  whole_array[i][j-1])
    return whole_array


def trace_shared_spliced(t1, t2, array):
    m = len(t1)
    n = len(t2)
    motif = ''
    while array[m][n] != 0:
        if array[m][n] == array[m][n-1]:
            n -= 1
        elif array[m][n] == array[m-1][n]:
            m -= 1
        else:
            motif += t1[m-1]
            m -= 1
            n -= 1
    motif = list(motif)
    motif.reverse()
    motif = ''.join(motif)
    return motif


def find_ssm(DNA1, DNA2):
    mainlst = []
    for _ in range(len(DNA1)+1):
        lst = []
        mainlst.append(lst)
        for _ in range(len(DNA2)+1):
            lst.append(0)
    for i in range(1, len(DNA1)+1):
        for j in range(1, len(DNA2)+1):
            if DNA1[i-1] == DNA2[j-1]:
                mainlst[i][j] = mainlst[i-1][j-1]+1
            else:
                mainlst[i][j] = max(mainlst[i-1][j], mainlst[i][j-1])
    m, n, lcmq = len(DNA1), len(DNA2), ''
    while mainlst[m][n] != 0:

        if mainlst[m][n] == mainlst[m-1][n]:
            m -= 1
        elif mainlst[m][n] == mainlst[m][n-1]:
            n -= 1
        else:
            lcmq += DNA1[m-1]
            m -= 1
            n -= 1
    lcmq = list(lcmq)
    lcmq.reverse()
    lcmq = ''.join(lcmq)
    return lcmq

data = reader('data_file')
t = trace_shared_spliced(data[0], data[1], shared_spliced_array(data[0], data[1]))
print t
print find_ssm(data[0], data[1])
if spliced_motif(data[0], t) and spliced_motif(data[1], t):
    print 'True'
else:
    print 'false'