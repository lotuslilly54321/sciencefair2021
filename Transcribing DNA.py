def transcribing_dna(string):
    new_string = ""
    for char in string:
        if char == 'T':
            new_string += 'U'
        else:
            new_string += char
    return new_string


print(transcribing_dna('CAAGTCGTTCACTGGTGCGTGCGTAGACTACGGGTTATCAGATGTAGCTAGAAGCACACTCGTGCACTCAGACAATCCGCCTCGATAGATTGAAGCATCGACTAATACCGGTATCAACCTAGAGTTCACGTAATCCAACCATCTAAGACGATCCCCCGGTTGTTCGAATCACGCTTTGTAACACATGTGATGACCGTTGAGGCCCCTTAGGCTGGGGTTTTGTGTAATTGTAAGACGAAGCCTGGCCAACTTACCATCAGCCTACCTCTTGGCTGATTAACATACGCTGAAAATTTAACGTCGTCCTAGTTGTCACGAGAAGGTGCTTTAAATAGTAGCTTTTACCATCCTCTACTCAGCGGGGTCACAGTTAACCCCGTAAAGGGGGTAATCACCAAACCGTACTAACCAAAAGTGCTAAGATACTTTGAGTGAGGGTTCTGTTCTTCATGAAGTCCAAACCTGGCGGCAGCGCGATTTGACTCCTTTTGGGCGGCATTAGCTGAGCAAAGCACCATCCACTCCGATTACTGGCGATAGGCTCCTTATATCCACCTGAACCCTCTCAGACACCGCTCTTCACTACATTACCCTGATTGGATGAGCTCGCTTTAAATGTTGCTGTTAATGGGGCGGGCCGTCGTCTATTCGCTTTGGGTCCAACGCAATTCAGGGTAATCGCACGGTGGTCTCCTGTAGCAGAATCCCTGTCCACATTGCTGCTCCTAAGGGAGTAACCCAGAGCAGGAAGCCTTACGTCCGGATTAAAGCTAGAGTTACAGCTGGGAGAGTTTCATTACGTAGGTACCAAGGGGGTCGGCGGCAACATGAACATGAATCTGCGGCATTACAGGACCATGTTAGGCGGGATTGATCACCATATAACTCTATTCATTTGGCCACTCGGT'))