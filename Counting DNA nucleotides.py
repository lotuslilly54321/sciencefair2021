def counting_dna_nucleotides(string):
    return count_char('A',string), count_char('C',string), count_char('G',string), count_char('T',string)


def count_char(char, string):
    count = 0
    for c in string:
        if c == char:
            count += 1
    return count


print counting_dna_nucleotides("")