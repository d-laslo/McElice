import numpy as np


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

def mul(x1, x2, pol):
    if x1 == 0 or x2 == 0:
            return 0
        
    bx1 = bin(x1)[2:][::-1]
    indexes = [i for i in range(len(bx1)) if bx1[i] == '1']
    res = 0
    for i in indexes:
        res ^= x2 << i

    return norm(res, pol)


def div(divided, divider):
    quotient = 0
    ldivider = len(bin(divider)[2:])
    while divided >= divider:
        ldivided = len(bin(divided)[2:])
        quotient ^= 1 << (ldivided - ldivider)
        divided ^= divider << (ldivided - ldivider)
    return quotient


def pow(base, exponenta, pol):
    itr = bin(exponenta)[2:][::-1]
    tmp = base
    result = 1
    for i in itr:
        if i == '1':
            result = mul(result, tmp, pol)
        tmp = mul(tmp, tmp, pol)
    return result


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

def get_indexes(x):
    x = bin(x)[2:][::-1]
    ind = []
    for i in range(len(x)):
        if x[i] == '1':
            ind.append(i)

    return ind


def transpose(matrix):
    return list(np.array(matrix).transpose())


def binary_transpose(matrix, _len):
    bmatrix = []
    for i in matrix:
        bmatrix.append(list('0' * (_len - len(bin(i)[2:])) + bin(i)[2:]))

    bmatix_T = list(np.array(bmatrix).transpose())
    matrix_T = []
    for i in bmatix_T:
        matrix_T.append(int(''.join(i), 2))

    return matrix_T


def sum(a):
    res = 0
    for i in a:
        res ^= i

    return res


# def mul_vectors(v1, v2, pol):
#     if len(v1) != len(v2):
#         raise Exception("Mul vectors")

#     return


def mul_matrix(m1, m2, pol):
    # m1 = [
    #     [1, 5, 12, 4, 10, 3, 1, 3, 4, 2, 2, 10], 
    #     [8, 11, 13, 9, 1, 0, 9, 3, 13, 14, 12, 11]
    # ]

    # m2 = [
    #     [0, 0, 1, 1],
    #     [1, 1, 1, 1],
    #     [1, 1, 0, 1],
    #     [0, 1, 1, 0],
    #     [1, 1, 1, 1],
    #     [0, 0, 0, 1],
    #     [1, 0, 0, 0],
    #     [0, 1, 0, 1],
    #     [0, 1, 0, 0],
    #     [1, 0, 0, 0],
    #     [0, 0, 0, 1],
    #     [0, 0, 1, 0]
    # ]
    # m1 = convert_to_binary_matrix(m1, 4)
    m2 = transpose(m2)

    result = []
    for m1_row in m1:
        result.append([])
        for m2_row in m2:
            result[-1].append(sum([mul(x1, x2, pol) for x1, x2 in zip(m1_row, m2_row)]))

    return result


def convert_to_binary_matrix(matrix, element_len):
    # matrix = [
    #     [1, 5, 12, 4, 10, 3, 1, 3, 4, 2, 2, 10], 
    #     [8, 11, 13, 9, 1, 0, 9, 3, 13, 14, 12, 11]
    # ]
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



if __name__ == '__main__':
    pass