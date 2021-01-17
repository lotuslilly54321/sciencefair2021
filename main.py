from utils import read_fasta_data, print_in_color
from config import comparisons


def cost_of_char_substitution(char_in_string_1, char_in_string_2):
    """
    Return the cost of substituting a char in the string 1 with the char in string 2
    :param char_in_string_1: character in string 1
    :param char_in_string_2: character in string 2
    :return: 0 if the characters are the same, else 1
    """
    if char_in_string_1 == char_in_string_2:
        return 0
    return 1


def edit_distance(string_1, string_2):
    """
    Algorithm that uses dynamic programming to recursively calculate the optimal number of substitutions,
    insertions or deletions to string_1 and string_2,required to match the two strings

    :param string_1:  first string
    :param string_2:  second string
    :return: count of the optimal number of substitutions
    """
    gap = 1
    edit_distance_matrix = []
    # create 2D array with a length and width one greater than each string
    for i in range(len(string_1) + 1):
        matrix_row = []
        for j in range(len(string_2) + 1):
            matrix_row.append(0)
        edit_distance_matrix.append(matrix_row)
    # initialize array for initial costs associated with starting one string at a location other than the 0th index
    for i in range(1, len(string_1) + 1):
        edit_distance_matrix[i][0] = i * gap
    for j in range(1, len(string_2) + 1):
        edit_distance_matrix[0][j] = j * gap

    # calculate the costs of insertion, deletion, and substitution and record the minimum value
    for i in range(1, len(string_1) + 1):
        for j in range(1, len(string_2) + 1):
            edit_distance_matrix[i][j] = min(
                cost_of_char_substitution(string_1[i - 1], string_2[j - 1]) + edit_distance_matrix[i - 1][j - 1],
                gap + edit_distance_matrix[i - 1][j],
                gap + edit_distance_matrix[i][j - 1])
    return edit_distance_matrix[len(string_1)][len(string_2)]


def edit_array(string_1, string_2):
    gap = 1
    edit_distance_matrix = []
    # create 2D array with a length and width one greater than each string
    for i in range(len(string_1) + 1):
        matrix_row = []
        for j in range(len(string_2) + 1):
            matrix_row.append(0)
        edit_distance_matrix.append(matrix_row)
    # initialize array for initial costs associated with starting one string at a location other than the 0th index
    for i in range(1, len(string_1) + 1):
        edit_distance_matrix[i][0] = i * gap
    for j in range(1, len(string_2) + 1):
        edit_distance_matrix[0][j] = j * gap

    # calculate the costs of insertion, deletion, and substitution and record the minimum value
    for i in range(1, len(string_1) + 1):
        for j in range(1, len(string_2) + 1):
            edit_distance_matrix[i][j] = min(
                cost_of_char_substitution(string_1[i - 1], string_2[j - 1]) + edit_distance_matrix[i - 1][j - 1],
                gap + edit_distance_matrix[i - 1][j],
                gap + edit_distance_matrix[i][j - 1])
    i = len(string_1)
    j = len(string_2)
    string_1_edited = ''
    string_2_edited = ''
    while i > 0 and j > 0:
        above = edit_distance_matrix[i - 1][j]
        left = edit_distance_matrix[i][j - 1]
        diagonal = edit_distance_matrix[i - 1][j - 1]
        if above < left and above < diagonal:
            string_2_edited += '-'
            string_1_edited += string_1[i - 1]
            i -= 1
        elif left < diagonal:
            string_2_edited += string_2[j - 1]
            string_1_edited += '-'
            j -= 1
        else:
            if string_1[i - 1] == string_2[j - 1]:
                string_1_edited += string_1[i - 1]
                string_2_edited += string_2[j - 1]
            else:
                string_1_edited += string_1[i - 1].lower()
                string_2_edited += string_2[j - 1].lower()
            i -= 1
            j -= 1
    string_1_edited_as_list = list(string_1_edited)
    string_2_edited_as_list = list(string_2_edited)
    string_1_edited_as_list.reverse()
    string_2_edited_as_list.reverse()
    new_string_1_edited = ''.join(string_1_edited_as_list)
    new_string_2_edited = ''.join(string_2_edited_as_list)
    return new_string_1_edited, new_string_2_edited


def main():
    """
    Compare the different FASTA files as stored in the comparisons config
    """
    for comparison in comparisons:
        desc = comparison["desc"]
        file_1 = comparison["file_1"]
        file_2 = comparison["file_2"]
        string_1 = read_fasta_data(file_1)
        string_2 = read_fasta_data(file_2)
        distance = edit_distance(string_1, string_2)
        print "----------------------------------------------------------------------------------"
        print desc
        print "----------------------------------------------------------------------------------"
        edited_string_1, edited_string_2 = edit_array(string_1, string_2)
        print_in_color(edited_string_1, 'green', 'red')
        print_in_color(edited_string_2, 'green', 'blue')
        length = len(edited_string_1)
        print "Percent Edit Distance={0}".format(float(distance)/length)
        print ""


if __name__ == "__main__":
    main()
