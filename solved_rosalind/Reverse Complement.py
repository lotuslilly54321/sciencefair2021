def reverse_complement(string):
    new_string = ""
    for char in string:
        if char == 'A':
            new_string = "T" + new_string
        if char == 'C':
            new_string = "G" + new_string
        if char == 'T':
            new_string = "A" + new_string
        if char == 'G':
            new_string = "C" + new_string
    return new_string


print reverse_complement('AGTTTTCTCCTGCACTCGAACTCCGGAATGGGCGCGAGTTTTCAAGCATTGTTATTGAAGGAACAAATGAAGAAGACCGTAGAGACCTACTAACCCGCTTAGAAGGTAAAGTTATGTGACACCTTGGGCTTTGTAGGTAAGGGTCTCACCGAAAGAGGGCGACTACAACTGTCCATAAGCGCTAAGCGTGAGCGCGAAATCGTATGATTCCTAGGATCCCGCGTATCTCGACACAGTAAAGCGCAACGGTCTTTCAGCGGGGGTGTCGCACAATCACAGTGTCCGAAAGGGTAGGTACAGTAGACGCTGATTTTCCCAAACTAATTCTACTTTTGGTTTCAGACAGGGTCATATCCATGACAAGCGTCTGGGTTTCTTGGAGAAGAATCCACGCTCCAGCGCCGCGGGAATCCAGGAGCACATACGGGTTGTGTACCCGTCGGAGAAACGAGTGCTAAGTTGTACTGTGTGACCGGAGTTATCACAATCGACGTCCATTGTTCTTACTTTACAGGAGCACTCAGCAGCCCGGTATTTAGACCACACCCCACTTTGCTCGCGCACCCTTAGTCGTGTGTCTACTAGTCAGATAGATACCGGACCACCCCCATGGACCCCCATTGAGAGCGCAACGCCCTTATAACTGAATATCACAGCTACCAAGAGCTTCTCGCCAGGTAAATTAACTCCAGGAGGGCGAATACGCCTTGAACCCTAGACTCCAGTAGAGGGCGGGGATAGTGGGTGCTCGCTATGGTCTATCTATGCGTACAGCGTCTGAGTTTCAGGATAAGAGGCAGGTTCCGGCAACTGCGCTAACGTACCAGGGCTAAATCACTATCGAGCTTGGGCCAGAGGGCCAGAGGGTCTTAGTTGTAACAAGACCAACCATGAAGTGGTTAATGCAACGGGTCGTTACTATG')