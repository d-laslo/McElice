import numpy as np
from copy import copy
import re


def read(path):
    with open(path, 'r') as file:
        text = file.read()
        return text

def write(path, text):
    with open(path, 'w') as file:
        file.write(text)


def norm(x, pol):
    bx = bin(x)[2:]
    bpol = bin(pol)[2:]

    while len(bx) >= len(bpol):
        shift = len(bx) - len(bpol)
        x ^= pol << shift
        bx = bin(x)[2:]
    return x


def mod(x, pol):
        bx = bin(x)[2:]
        bpol = bin(pol)[2:]

        while x >= pol:
            shift = len(bx) - len(bpol)
            x ^= pol << shift
            bx = bin(x)[2:]
        return x


# множення елементів із поля яке генерує поліном pol
# pol заданий як десяткове число
def mul(x1, x2, pol):
    if x1 == 0 or x2 == 0:
            return 0
        
    bx1 = bin(x1)[2:][::-1]
    indexes = [i for i in range(len(bx1)) if bx1[i] == '1']
    res = 0
    for i in indexes:
        res ^= x2 << i

    return norm(res, pol)


# множення поліномів коефіцієнти якого є елементами із поля яке генерує поліном pol
# pol заданий як десяткове число
def mul_p(multiplier1, multiplier2, pol):
    if len(multiplier1) > len(multiplier2):
        multiplier1, multiplier2 = multiplier2, multiplier1

    mult1_s = len(multiplier1)
    mult2_s = len(multiplier2)
    res = [0 for i in range(mult1_s + mult2_s - 1)]

    for i in range(mult2_s):
        tmp = copy(multiplier1)
        for j in range(len(tmp)):
            tmp[j] = mul(tmp[j], multiplier2[i], pol)

        tmp = tmp + [0 for i in range(mult1_s + mult2_s - i - 1 - len(tmp))]
        tmp = tmp[::-1]

        for j in range(len(tmp)):
            res[j] ^= tmp[j]
    res = res[::-1]
    return res


# ділення двох елементів із поля, що генерує поліном pol
# pol заданий як десяткове число
def div(divided, divider):
    quotient = 0
    ldivider = len(bin(divider)[2:])
    while divided >= divider:
        ldivided = len(bin(divided)[2:])
        quotient ^= 1 << (ldivided - ldivider)
        divided ^= divider << (ldivided - ldivider)
    return quotient


# ділення поліномів коефіцієнти якого є елементами із поля яке генерує поліном pol
# pol заданий як десяткове число
def div_p(divided, divider, pol):
    if len(divider) > len(divided):
        return [[0], divided]

    quotient = [0 for i in range(len(divided))]
    divider_l = len(divider)

    inv = mod_inverse(divider[0], pol)
    for i in range(divider_l):
        divider[i] = mul(divider[i], inv, pol)


    while len(divided) >= divider_l:
        divided_l = len(divided)

        mult = divided[0]
        quotient[divided_l - divider_l] = mul(mult, inv, pol)

        tmp = copy(divider)
        for i in range(len(tmp)):
            tmp[i] = mul(tmp[i], mult, pol)

        tmp = tmp + [0 for i in range(divided_l - divider_l)]
        for i in range(divided_l):
            divided[i] ^= tmp[i]

        while len(divided) > 0 and divided[0] == 0:
            divided.pop(0)
    
    

    quotient = quotient[::-1]
    while len(quotient) > 0 and quotient[0] == 0:
        quotient.pop(0)
    if len(divided) == 0:
        divided = [0]
    return [quotient, divided]


# піднесення до степення елементів поля, що генерує поліном pol
# pol заданий як десяткове число
def pow(base, exponenta, pol):
    itr = bin(exponenta)[2:][::-1]
    tmp = base
    result = 1
    for i in itr:
        if i == '1':
            result = mul(result, tmp, pol)
        tmp = mul(tmp, tmp, pol)
    return result


# обернений за модулем елемент із поля яке генерує поліном pol
# pol заданий як десяткове число
def mod_inverse(value, module):
    v = [0, 1]
    rem = 1
    a = module

    while rem > 0:
        tmp = div(a, value)
        v[0] ^= mul(v[1], tmp, module)
        v[0], v[1] = v[1], v[0]
        rem = mod(a, value)
        a = value
        value = rem
    return v[0]


# обернений за поліномом g поліном, елементи якого із поля яке генерує поліном pol
# pol заданий як десяткове число
# g заданий як list; g[0] -- найстарший елемент
def mod_inverse_p(element, g, pol):
    if len(element) == 1:
        return [mod_inverse(element[0], pol)]

    v = [[0], [1]]
    rem = None
    value = element
    a = copy(g)


    while rem == None or len(rem) > 1:
        t = div_p(copy(a), copy(value), pol)
        tmp = t[0]
        rem = t[1]

        tmp_mul = mul_p(copy(v[1]), copy(tmp), pol)
        v[0] = v[0] + [0 for i in range(len(tmp_mul) - len(v[0]))]
        for i in range(len(v[0])):
            v[0][i] ^= tmp_mul[i]
        v[0], v[1] = v[1], v[0]
        a = copy(value)
        value = copy(rem)
    
    inv = mod_inverse(rem[0], pol)
    for i in range(len(v[1])):
        v[1][i] = mul(inv, v[1][i], pol)
    return v[1]


# повертає індекси не нулбових елементів у двійковому записі числа x
def get_indexes(x):
    x = bin(x)[2:][::-1]
    ind = []
    for i in range(len(x)):
        if x[i] == '1':
            ind.append(i)

    return ind


# транспонує матрицю
def transpose(matrix):
    return list(np.array(matrix).transpose())


# транспонує біти матриці записаної у десятковому записі
# _len вказує максимальну довжину елемента у бітовому записі
def binary_transpose(matrix, _len):
    bmatrix = []
    for i in matrix:
        bmatrix.append(list('0' * (_len - len(bin(i)[2:])) + bin(i)[2:]))

    bmatix_T = list(np.array(bmatrix).transpose())
    matrix_T = []
    for i in bmatix_T:
        matrix_T.append(int(''.join(i), 2))

    return matrix_T


# сума елементів із поля характеристики 2
def sum(a):
    res = 0
    for i in a:
        res ^= i

    return res


def convert_to_binary_matrix(matrix, element_len):
    matrix = transpose(matrix)
    bmatrix = []
    for row in matrix:
        bmatrix.append([int(em) for em in list(''.join(['0' * (element_len - len(bin(i)[2:])) + bin(i)[2:] for i in row]))])
    bmatrix = transpose(bmatrix)

    return bmatrix


def index_max_len(elements):
    ind = 0
    max_len = 0
    for i in range(len(elements)):
        if max_len < len(bin(elements[i])):
            max_len = len(bin(elements[i]))
            ind = i
    return ind


def convert2num(pol):
    pows = [int(i) for i in re.findall(r'\w\^(\d*)', pol)]
    res = 1
    for i in pows:
        res ^= (1 << i)

    return res

def convert2vector(pol):
    pows = [int(i) for i in re.findall(r'\w\^(\d*)', pol)]
    res = [0 for i in range(max(pows) + 1)]
    for i in pows:
        res[i] = 1
    res[0] = 1
    return res[::-1]
 

if __name__ == '__main__':
    pass
    print(convert2vector('x^5+x^3+1'))