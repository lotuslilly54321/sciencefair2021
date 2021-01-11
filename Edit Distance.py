from SeqIO import read_fasta_file_modif


def reader(string):
    data = read_fasta_file_modif(string)
    return data

#
# def shared_spliced_array(t1, t2):
#     whole_array = []
#     for i in range(len(t1)+1):
#         lst = []
#         for j in range(len(t2)+1):
#             lst.append(0)
#         whole_array.append(lst)
#     for i in range(1, len(t1)+1):
#         for j in range(1, len(t2)+1):
#             if t1[i-1] == t2[j-1]:
#                 whole_array[i][j] = whole_array[i-1][j-1] + 1
#             else:
#                 whole_array[i][j] = max(whole_array[i-1][j],  whole_array[i][j-1])
#     return whole_array
#
#
# def trace_shared_spliced(t1, t2, array):
#     m = len(t1)
#     n = len(t2)
#     motif = ''
#     while array[m][n] != 0:
#         if array[m][n] == array[m][n-1]:
#             n -= 1
#         elif array[m][n] == array[m-1][n]:
#             m -= 1
#         else:
#             motif += t1[m-1]
#             m -= 1
#             n -= 1
#     motif = list(motif)
#     motif.reverse()
#     motif = ''.join(motif)
#     return motif


def mms(short, long):
    mms = []
    for i in range(len(long)-len(short)+1):
        mm = 0
        for j in range(len(short)):
            if long[i+j] == short[j]:
                pass
            else:
                mm += 1
        mms.append(mm)
    return mms


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


# def edit_array(s1, s2):
#     gap = 1
#     A = []
#     B = []
#     C = []
#     # create 2D array
#     for i in range(len(s1)+1):
#         a = []
#         b = []
#         c = []
#         for j in range(len(s2)+1):
#             a.append(0)
#             b.append((0, 0))
#             c.append('')
#         A.append(a)
#         B.append(b)
#         C.append(c)
#
#     # initialize array for base cases
#     for i in range(1, len(s1)+1):
#         A[i][0] = i * gap
#     for j in range(1, len(s2)+1):
#         A[0][j] = j * gap
#
#     # calculate the costs
#     for i in range(1, len(s1)+1):
#         for j in range(1, len(s2)+1):
#             minimum = min(
#                             cost(s1[i-1], s2[j-1]) + A[i-1][j-1],
#                             gap + A[i-1][j],
#                             gap + A[i][j-1])
#             A[i][j] = minimum
#             if minimum == gap + A[i][j-1]:
#                 B[i][j] = (i, j-1, 'a')
#             elif minimum == gap + A[i-1][j]:
#                 B[i][j] = (i-1, j, 'b')
#             elif minimum == cost(s1[i - 1], s2[j - 1]) + A[i - 1][j - 1]:
#                 B[i][j] = (i - 1, j - 1, 'p')
#     previous = B[i][j]
#     edited_a = list(str1)
#     edited_b = list(str2)
#     while True:
#         prev_i, prev_j, edit = previous
#         if prev_i == 0 and prev_j == 0:
#             break
#         if edit == 'a':
#             edited_a.insert(prev_i-1, '-')
#         if edit == 'b':
#             edited_b.insert(prev_j-1, '-')
#         previous = B[prev_i][prev_j]
#     print ''.join(edited_a)
#     print ''.join(edited_b)
#     return edited_a, edited_b


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
            a_edited += str1[i-1]
            i -= 1
        elif left < diagonal:
            b_edited += str2[j-1]
            a_edited += '-'
            j -= 1
        else:
            a_edited += str1[i-1]
            b_edited += str2[j-1]
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





data = reader('data_file')
str1 = data[0]
str2 = data[1]


# print mms('AAA', 'GGGAAGAAAG')
# print edit_distance('AAA', '0123AAA7890')
print edit_distance(data[0], data[1])
print edit_array(data[0], data[1])