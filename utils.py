from __future__ import print_function
from termcolor import colored

from Bio import SeqIO
import sys
import json


#
# Reads fasta file data
#
def read_fasta_data(string):
    data_arr = read_fasta_file_as_data_arr(string)
    fasta_data = ''
    for data in data_arr:
        fasta_data += data
    return fasta_data

#
#
#
def read_fasta_file(filename):
    data = ""
    with open(filename, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            # identifier = record.id
            # description = record.description
            data = data + record.seq
    data = data.strip('Z')
    return data


def read_fasta_file_as_data_arr(filename):
    data_set = []
    with open(filename, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            # identifier = record.id
            # description = record.description
            data = str(record.seq)
            data = data.strip('Z')
            data_set.append(data)
    return data_set
#
#
#
def naive_approximate(gene, chromosome, max_hamming_distance):
    occurs = []
    for i in xrange(0, len(chromosome) - len(gene) + 1):  # for	all	alignments
        nmm = 0
        for j in xrange(0, len(gene)):
            if chromosome[i + j] != gene[j]:  # does	it	match?
                nmm += 1
                if nmm > max_hamming_distance:
                    break  # exceeded	maximum	distance
        if nmm <= max_hamming_distance:
            occurs.append((i, nmm))
    return occurs


def print_in_color(s, standard_color, point_mutation_color):
    colored_s = ''
    for c in s:
        if c.islower():
            colored_s += colored(c.upper(), point_mutation_color)
        else:
            colored_s += colored(c, standard_color)
    print(colored_s)


class TerminalColors:
    def __init__(self):
        pass
    HEADER = '\033[95m'
    OK_BLUE = '\033[94m'
    OK_CYAN = '\033[96m'
    OK_GREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    END_COLOR = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#
# chr1 = read_fasta_file('data_file/chr1.FASTA')
# chr2 = read_fasta_file('data_file/chr2.FASTA')
# chr3 = read_fasta_file('data_file/chr3.FASTA')
# chr4 = read_fasta_file('data_file/chr4.FASTA')
# rubens = read_fasta_file('data_file/UP000000724.fasta')
#
# chromosomes_dictionary = {
#     "chromosome1": chr1,
#     "chromosome2": chr2,
#     "chromosome3": chr3,
#     "chromosome4": chr4,
# #    "rubens": rubens,
# }
#
# pcbAB = read_fasta_file('data_file/A6N339.fasta')
# pcbC = read_fasta_file('data_file/P08703.fasta')
# penDE = read_fasta_file('data_file/P15802.fasta')
#
# genes_dictionary = {
#     "pcbAB": pcbAB,
#     "pcbC": pcbC,
#     "penDE": penDE,
# }
#
#
# def protein_parts(protein):
#     parts = []
#     for i in range(0, len(protein), 10):
#         start_index = i
#         end_index = min(i + 10, len(protein))
#         parts.append(str(protein[start_index:end_index]))
#     return parts
#
#
# # attempt to find using hamming distance on the entire protein
#
# max_hamming_distance = 3
#
# genes_found_dictionary = {}
#
# details_file = open('PenicilliumSearchDetails.txt', 'w')
#
# for chromosome_key, chromosome_value in chromosomes_dictionary.items():
#     for gene_key, gene_value in genes_dictionary.items():
#         occurrences = naive_approximate(gene_value, chromosome_value, max_hamming_distance)
#         for occurrence in occurrences:
#             location, num_differences = occurrence
#             genes_found_dictionary[gene_key] = location
#             print('found gene {} in chromosome {} at location {} with {} differences'.format(
#                 gene_key, chromosome_key, location, num_differences))
#             print('found gene {} in chromosome {} at location {} with {} differences'.format(
#                 gene_key, chromosome_key, location, num_differences), file=details_file)
#
# # for any protein not matched, its possible that there was a shift (insertion or deletion)
# # that caused the hamming distance algorithm to fail
#
# genes_not_found_list = []
# for gene_key, gene_value in genes_dictionary.items():
#     if gene_key not in genes_found_dictionary:
#         genes_not_found_list.append(gene_key)
#
# print('Genes not found {}'.format(genes_not_found_list))
# print('Genes not found {}'.format(genes_not_found_list), file=details_file)
#
# print('Attempting to run hamming distance on shorter substrings')
# print('Attempting to run hamming distance on shorter substrings', file=details_file)
#
# max_allowed_parts_mismatch_percent = 2
#
# # attempt to run hamming distance on a shorter substrings
# for chromosome_key, chromosome_value in chromosomes_dictionary.items():
#     for gene_key, gene_value in genes_dictionary.items():
#         overall_mismatch = 0
#         if gene_key in genes_not_found_list:
#             parts_matched = []
#             parts_not_matched = []
#             gene_part_values = protein_parts(gene_value)
#             max_allowed_parts_mismatch = int(len(gene_part_values) * 1.0 * max_allowed_parts_mismatch_percent) / 100.0
#             parts_mismatch_count = 0
#             search_start_index = 0
#             min_location = sys.maxint
#             for gene_part_value in gene_part_values:
#                 occurrences = naive_approximate(gene_part_value, chromosome_value[search_start_index:],
#                                                 max_hamming_distance)
#                 # print('occurrences={}'.format(occurrences))
#                 min_occurrence_location = sys.maxint
#                 for occurrence in occurrences:
#                     location, num_differences = occurrence
#                     location = location + search_start_index
#                     min_occurrence_location = min(min_occurrence_location, location)
#                     overall_mismatch += num_differences
#                     # print('**debug** found part_value {} in chromosome {} at location {}'.format(gene_part_value,
#                     #                                                                              chromosome_key,
#                     #                                                                              location))
#                 search_start_index = max(search_start_index, min_occurrence_location)
#                 min_location = min(min_location, min_occurrence_location)
#                 if len(occurrences) == 0:
#                     # print('**debug** part_value {} not found in chromosome {}'.format(
#                     # gene_part_value,
#                     # chromosome_key))
#                     parts_mismatch_count += 1
#                     parts_not_matched.append(gene_part_value)
#                 else:
#                     parts_matched.append(gene_part_value)
#             if parts_mismatch_count <= max_allowed_parts_mismatch:
#                 mismatch_percent = ((100.0 * overall_mismatch)+(len(parts_not_matched)*10))/len(gene_value)
#                 print('found gene {} in chromosome {} at location {}'.format(
#                     gene_key, chromosome_key, min_location, ))
#                 print('found gene {} in chromosome {} at location {}'.format(
#                     gene_key, chromosome_key, min_location), file=details_file)
#                 print('***SUMMARY**')
#                 print('parts_matched count={} parts_not_matched count={}'.format(
#                     len(parts_matched), len(parts_not_matched)))
#                 print('parts_matched count={} parts_not_matched count={}'.format(
#                     len(parts_matched), len(parts_not_matched)), file=details_file)
#                 print('***DETAILS are in PenicilliumSearchDetails.txt**')
#                 print('PARTS NOT MATCHED', file=details_file)
#                 print(json.dumps(parts_not_matched), file=details_file)
#                 print('PARTS MATCHED', file=details_file)
#                 print(json.dumps(parts_matched), file=details_file)
#                 print('MISMATCH PERCENT {}'.format(mismatch_percent))
#                 print('MISMATCH PERCENT {}'.format(mismatch_percent), file=details_file)
#
# details_file.close()
