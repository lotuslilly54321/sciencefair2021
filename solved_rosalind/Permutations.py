from itertools import permutations
import math


def starter(num):
    starter_string = ''
    for i in range(int(num)):
        starter_string = starter_string + str(i + 1)
    return starter_string

def all_permutations(num):
    starter_string = starter(num)
    print(math.factorial(num))
    for p in permutations(starter_string):
        output = ""
        for item in p:
            output += ' ' + item
        print output[1:]


all_permutations(7)