from SeqIO import read_fasta_file_modif


def reader(string):
    data = read_fasta_file_modif(string)
    return data


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
        for j in range(1,len(s2)+1):
            A[i][j] = min(
                            cost(s1[i-1], s2[j-1]) + A[i-1][j-1],
                            gap + A[i-1][j],
                            gap + A[i][j-1])
    return A[len(s1)-1][len(s2)-1]


def edit_array(s1, s2):
    gap = 1
    A = []
    B = []
    # create 2D array
    for i in range(len(s1)+1):
        a = []
        b = []
        for j in range(len(s2)+1):
            a.append(0)
            b.append((0, 0))
        A.append(a)
        B.append(b)

    # initialize array for base cases
    for i in range(1, len(s1)+1):
        A[i][0] = i * gap
    for j in range(1, len(s2)+1):
        A[0][j] = j * gap

    # calculate the costs
    for i in range(1, len(s1)+1):
        for j in range(1,len(s2)+1):
            minimum = min(
                            cost(s1[i-1], s2[j-1]) + A[i-1][j-1],
                            gap + A[i-1][j],
                            gap + A[i][j-1])
            A[i][j] = minimum
            if minimum == cost(s1[i-1], s2[j-1]) + A[i-1][j-1]:
                B[i][j] = (i-1, j-1)
            elif minimum == gap + A[i-1][j]:
                B[i][j] = (i-1, j)
            elif minimum == gap + A[i][j-1]:
                B[i][j] = (i, j-1)
    paths = []
    paths.append((i,j))
    previous = B[i][j]
    while True:
        prev_i, prev_j = previous
        if prev_i == 0 and prev_j == 0:
            break
        paths.append(previous)
        previous = B[prev_i][prev_j]
    return paths


def edit_alignment(str1, str2):
    path = edit_array(str1, str2)
    built_strings = ['', '']
    i = len(str1) - 1
    j = len(str2) - 1
    while len(path) > 0:
        if path[len(path)-1] == 'p':
            built_strings[0] = str1[i] + built_strings[0]
            built_strings[1] = str2[j] + built_strings[1]
            i -= 1
            j -= 1
        elif path[len(path)-1] == 's':
            built_strings[0] = '-' + built_strings[0]
            built_strings[1] = str2[j] + built_strings[1]
            j -= 1
        elif path[len(path)-1] == 'd':
            built_strings[0] = str1[i] + built_strings[0]
            built_strings[1] = '-' + built_strings[1]
            i -= 1
        print path
        path = path[1:]
        print built_strings




data = reader('data_file')
str1 = data[0]
str2 = data[1]


# print mms('AAA', 'GGGAAGAAAG')
# print edit_distance('AAA', '0123AAA7890')

print edit_alignment(data[0], data[1])