from SeqIO import read_fasta_file_modif


def reader(string):
    data = read_fasta_file_modif(string)
    return data


def cost(c_i, c_j):
    if c_i == c_j:
        return 0
    return 1


def naive_mm(str1, str2):
    mm = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            mm += 1
    return mm


def all_indexes(char, list):
    indexes = []
    for i in range(len(list)):
        if list[i] == char:
            indexes.append(i)
    return indexes


def edit_distance(s1, s2):
    gap = 1
    A = []
    # create 2D array
    for i in range(len(s1)+1):
        a = []
        for j in range(len(s2)+1):
            a.append(0)
        A.append(a)
    # initialize array for base cases
    for i in range(1, len(s1)+1):
        A[i][0] = i * gap
    for j in range(1, len(s2)+1):
        A[0][j] = j * gap

    # calculate the costs
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            A[i][j] = min(
                            cost(s1[i-1], s2[j-1]) + A[i-1][j-1],
                            gap + A[i-1][j],
                            gap + A[i][j-1])
    return A[len(s1)][len(s2)]


def edit_array(s1, s2):
    gap = 1
    A = []
    # create 2D array
    for i in range(len(s1)+1):
        a = []
        for j in range(len(s2)+1):
            a.append(0)
        A.append(a)
    # initialize array for base cases
    for i in range(1, len(s1)+1):
        A[i][0] = i * gap
    for j in range(1, len(s2)+1):
        A[0][j] = j * gap

    # calculate the costs
    for i in range(1, len(s1)+1):
        for j in range(1, len(s2)+1):
            A[i][j] = min(
                            cost(s1[i-1], s2[j-1]) + A[i-1][j-1],
                            gap + A[i-1][j],
                            gap + A[i][j-1])
    i = len(s1)
    j = len(s2)
    a_edited = ''
    b_edited = ''
    while i > 0 and j > 0:
        above = A[i-1][j]
        left = A[i][j-1]
        diagonal = A[i-1][j-1]
        if above < left and above < diagonal:
            b_edited += '-'
            a_edited += s1[i-1]
            i -= 1
        elif left < diagonal:
            b_edited += s2[j-1]
            a_edited += '-'
            j -= 1
        else:
            if s1[i-1] == s2[j-1]:
                a_edited += s1[i-1]
                b_edited += s2[j-1]
            else:
                a_edited += s1[i - 1].lower()
                b_edited += s2[j - 1].lower()
            i -= 1
            j -= 1
    a_list = list(a_edited)
    b_list = list(b_edited)
    a_list.reverse()
    b_list.reverse()
    new_a = ''.join(a_list)
    new_b = ''.join(b_list)
    print new_a
    print new_b
    return new_a, new_b, naive_mm(new_a, new_b)


def edit_alignment(str1, str2):
    path = edit_array(str1, str2)
    built_strings = ['', '']
    i = len(str1) - 1
    j = len(str2) - 1





datas1 = reader('Hep E Burma.fasta')
data1 = ''
for data in datas1:
    data1 += data
string_1 = data1
datas2 = reader('Hep E Pakistan.fasta')
data2 = ''
for data in datas2:
    data2 += data
string_2 = data2


distance =  edit_distance(string_1, string_2)
print distance
length = len(edit_array(string_1, string_2)[0])
print distance, length, (float(distance)/length)*100
