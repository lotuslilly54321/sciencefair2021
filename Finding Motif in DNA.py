def finding_motif_in_dna(dna, motif):
    indexes = ''
    k = len(motif)
    for i in range(len(dna)-len(motif)+1):
        if dna[i:i+k] == motif:
            indexes += str(i+1) + ' '
    return indexes


print(finding_motif_in_dna('', ''))