from SeqIO import read_fasta_file_modif


def reader(string):
    data = read_fasta_file_modif(string)
    return data


codon_chart = {'UUU': 'F',
               'CUU': 'L',
               'AUU': 'I',
               'GUU': 'V',
               'UUC': 'F',
               'CUC': 'L',
               'AUC': 'I',
               'GUC': 'V',
               'UUA': 'L',
               'CUA': 'L',
               'AUA': 'I',
               'GUA': 'V',
               'UUG': 'L',
               'CUG': 'L',
               'AUG': 'M',
               'GUG': 'V',
               'UCU': 'S',
               'CCU': 'P',
               'ACU': 'T',
               'GCU': 'A',
               'UCC': 'S',
               'CCC': 'P',
               'ACC': 'T',
               'GCC': 'A',
               'UCA': 'S',
               'CCA': 'P',
               'ACA': 'T',
               'GCA': 'A',
               'UCG': 'S',
               'CCG': 'P',
               'ACG': 'T',
               'GCG': 'A',
               'UAU': 'Y',
               'CAU': 'H',
               'AAU': 'N',
               'GAU': 'D',
               'UAC': 'Y',
               'CAC': 'H',
               'AAC': 'N',
               'GAC': 'D',
               'UAA': 'Stop',
               'CAA': 'Q',
               'AAA': 'K',
               'GAA': 'E',
               'UAG': 'Stop',
               'CAG': 'Q',
               'AAG': 'K',
               'GAG': 'E',
               'UGU': 'C',
               'CGU': 'R',
               'AGU': 'S',
               'GGU': 'G',
               'UGC': 'C',
               'CGC': 'R',
               'AGC': 'S',
               'GGC': 'G',
               'UGA': 'Stop',
               'CGA': 'R',
               'AGA': 'R',
               'GGA': 'G',
               'UGG': 'W',
               'CGG': 'R',
               'AGG': 'R',
               'GGG': 'G'}


def transcribing_dna(string):
    new_string = ""
    for char in string:
        if char == 'T':
            new_string += 'U'
        else:
            new_string += char
    return new_string


def translate(rna):
    codons = []
    proteins = []
    protein_chain = ''
    for i in range(int(len(rna)/3)):
        codons.append(rna[i*3:(i*3)+3])
    for codon in codons:
        if codon_chart[codon] != 'Stop':
            proteins.append(codon_chart[codon])
        elif codon_chart[codon] == 'Stop':
            break
    for protein in proteins:
        protein_chain += protein
    return protein_chain


def rna_splicing(rna_strings):
    main = rna_strings[0]
    introns = rna_strings[1:]
    indexes = [0]
    exons = []
    final = main
    for intron in introns:
        if intron != 'Stop':
            final = final.replace(intron, '')
        else:
            break
    # for intron in introns:
    #     indexes.append(main.index(intron))
    #     indexes.append(main.index(intron)+len(intron))
    # for i in range(len(indexes)/2):
    #     exons.append(main[indexes[i*2]:indexes[(i*2)+1]])
    # exons.append(main[indexes[len(indexes)-1]:])
    # for exon in exons:
    #     final += exon
    return final




print(rna_splicing(reader('data_file')))


print translate(transcribing_dna(rna_splicing(reader('data_file'))))