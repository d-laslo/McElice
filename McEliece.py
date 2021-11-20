from GoppaCode import GoppaCode as GC
from source import *
from random import sample
import numpy as np
from copy import copy

class McEliece:
    def __init__(self, p, g, n = 0):
        self.__GC_flag = False
        if n != 0:
            self.__GC_flag = True
        self.__GC = GC(p, g, n)
        self.__G = []
        self.__H = []

        self.__S = []
        # self.__inv_S = []
        self.__P = [] #self.__gen_P()
        # self.__inv_P = [] #self.__calc_inv_P()
        self.__private_key = []
        self.__public_key = []


    @property
    def private_key(self):
        if len(self.__private_key) == 0:
            self.__private_key = self.__calc_private_key()
        return self.__private_key


    @property
    def publick_key(self):
        if len(self.__public_key) == 0:
            self.__public_key = self.__calc_public_key()
        return self.__public_key


    def __mul(self, x1, x2):
        return mul(x1, x2, self.__GC.p)


    def __dot(self, x1, x2):
        return np.dot(x1, x2).astype(np.uint64) % 2


    def __inverse(self, x):
        return mod_inverse(x, self.__GC.p)


    def __gen_P(self):
        sm = sample(range(self.__GC.n), self.__GC.n)
        P = np.zeros(self.__GC.n ** 2, np.uint64).reshape((self.__GC.n, self.__GC.n))
        for i in range(self.__GC.n):
            P[i, sm[i]] = 1
        return P


    def __calc_inv_P(self, P):
        if len(P) == 0:
            raise Exception()
        return np.array(P).T.astype(np.uint64)


    def __gen_S(self):
        t = 0
        S = []
        while t == 0:
            S = np.array(np.random.randint(0, 2, self.__GC.k ** 2))
            S = S.reshape((self.__GC.k, self.__GC.k))
            # S = np.array(np.random.randint(0, 2, 10000 ** 2))
            # S = S.reshape((10000, 10000))
            t = np.linalg.det(S) % 2
        return S


    def __calc_inv_S(self, S):
        if len(S) == 0:
            raise Exception()
        return np.linalg.inv(S).astype(np.uint64) % 2


    def __gen_key_parameters(self):
        if self.__GC_flag == None:
            raise Exception
        self.__S = self.__gen_S().astype(np.uint64)
        self.__P = self.__gen_P().astype(np.uint64)
        self.__G = np.array(self.__GC.G).astype(np.uint64)
        self.__H = np.array(self.__GC.H).astype(np.uint64)


    def __calc_private_key(self):
        if len(self.__S) == 0 or len(self.__P) == 0 or len(self.__G) == 0 or len(self.__H) == 0:
            self.__gen_key_parameters()
        return [self.__S, self.__G, self.__P, self.__H]


    def __calc_public_key(self):
        if len(self.__S) == 0 or len(self.__P) == 0 or len(self.__G) == 0:
            self.__gen_key_parameters()
        return [self.__dot(self.__dot(self.__S, self.__G), self.__P), self.__GC.t]


    def __convert_to_binary_vector(self, x, k):
        bx = '0' * (k - len(bin(x)[2:]))+ bin(x)[2:]
        if len(bx) > k:
            raise Exception()
        
        return [0 if i == '0' else 1 for i in list(bx)]


    def __convert_to_num(self, x):
        return int(''.join([str(i) for i in x]), 2)


    def __solution(self, equations):
        num_equtions = len(equations)
        row_len = len(equations[0])

        def mul_vector(vector, element):
            return [self.__mul(v, element) for v in vector]

        def length(vector):
            counter = 0
            for i in vector:
                if i != 0:
                    break
                counter += 1
            return len(vector) - counter

        def index(matrix, _len):
            for i in range(len(matrix)):
                if length(matrix[i]) == _len:
                    return i
            return -1

            
        for i in range(num_equtions):
            # max_element = max(equations[i:])
            ind = index(equations[i:], num_equtions + 1) + i
            max_element = equations[ind]
            len_max_element = num_equtions + 1 - i

            equations[ind] = equations[i]
            equations[i] = max_element
            
            norm_v = mul_vector(max_element, self.__inverse(max_element[0]))
            for j in range(i + 1, len(equations)):
                _len = length(equations[j])
                if len_max_element == _len:
                    t = np.array(mul_vector(norm_v, equations[j][row_len - _len]))
                    equations[j] = list(np.array(equations[j]) ^ np.array(mul_vector(norm_v, equations[j][row_len - _len])))

            equations = equations[::-1]


            x = []
            for i in range(num_equtions):
                ind = num_equtions - i - 1
                row = equations[i]
                row = mul_vector(row, self.__inverse(row[ind]))
                for j in range(ind + 1):
                    row.pop(0)

                b = row[-1]
                for a , c in zip(x[::-1], row):
                    b ^= self.__mul(a, c)
                x.append(b)
            return x


    def __get_mS(self, GG):
        equations = []

        def length(x):
            try:
                return self.__GC.n - list(x).index(1)
            except:
                return 0

        def find_indx(x, _len):
            for i in range(len(x)):
                if length(x[i]) == _len:
                    return i
            return -1

        for i in range(self.__GC.n):
            ind = find_indx(GG, self.__GC.n - i)
            max_element = copy(GG[ind])
            len_max_element = self.__GC.n - i
            
            GG[ind] = copy(GG[i])
            GG[i] = copy(max_element)

            for j in range(i + 1, self.__GC.n):
                if len_max_element == length(GG[j]):
                    GG[j] ^= max_element

        return [i[-1] for i in GG[:self.__GC.k]]
        # for row in GG:
        #     equations.append(int(''.join([str(em) for em in row]), 2))

        # # приводимо систему рівнянь до трапецеїдального вигляду
        # for i in range(len(equations)):
        #     # max_element = max(equations[i:])
        #     ind = index_max_len(equations[i:]) + i
        #     max_element = equations[ind]
        #     len_max_element = len(bin(max_element)) - 2

        #     equations[equations.index(max_element)] = equations[i]
        #     equations[i] = max_element

        #     for j in range(i + 1, len(equations)):
        #         if equations[j] > 0 and len_max_element == (len(bin(equations[j])) - 2):
        #             equations[j] ^= max_element

        # return [1 if bin(equations[i])[-1] == '1' else 0 for i in range(self.__GC.k)]
                
    
    def __syndrome(self, msg, H):
        H_T = transpose(convert_to_binary_matrix(H, self.__GC.m))

        # msg = np.array([0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1]).astype(np.uint64)
        tt = self.__dot(msg, H_T)

        def split(x, len_element):
            num_elements = len(x) // len_element
            return [
                int(''.join([
                    str(j) for j in x[i * len_element: (i + 1) * len_element]
                ]), 2) 
                for i in range(num_elements)
            ] 

        el = split(tt, self.__GC.m)
        if max(el) == 0:
            return 0
        
        matrix = []
        for i in range(self.__GC.t):
            matrix.append([el[j + i] for j in range(self.__GC.t + 1)])

        def sigma(_x, pol):
            x = 1
            sum = 0
            for i in pol:
                sum ^= self.__mul(x, i)
                x = self.__mul(x, _x)
            return sum

        pol = [1] + self.__solution(matrix)
        roots = [i for i in self.__GC.L if sigma(i, pol) == 0]

        inv_roots = [self.__inverse(r) for r in roots]
        indx = [self.__GC.L.index(el) for el in inv_roots]
        res = np.zeros(self.__GC.n).astype(np.uint64)
        for i in indx:
            res[i] = 1
        return res


    def encrypt(self, msg, public_key):
        Gh, t = public_key
        k = len(Gh)
        n = len(Gh[0])

        err = np.zeros(n, np.uint64)
        for i in sample(range(n), t):
            err[i] = 1

        bmsg = self.__convert_to_binary_vector(msg, k)
        # t = self.__dot(bmsg, self.__S)
        err = self.__dot(err, self.__P)
        return self.__convert_to_num(self.__dot(bmsg, Gh) ^ err)


    def decrypt(self, msg, private_key):
        S, G, P, H = private_key
        
        # домножуємо на оберену матрицю P
        bmsg = self.__convert_to_binary_vector(msg, self.__GC.n)
        bmsg = self.__dot(bmsg, self.__calc_inv_P(P))

        # знаходимо помилки
        err = self.__syndrome(bmsg, H)
        # # прибираємо помилки
        bmsg = bmsg ^ err

        # bmsg = bmsg[self.__GC.r:]
        GG = np.concatenate((transpose(G), bmsg.reshape(self.__GC.n, 1)), axis=1)
        bmsg = self.__get_mS(GG)
        # домножуємо на оберену матрицю P
        d_msg = self.__dot(bmsg, self.__calc_inv_S(S))

        return self.__convert_to_num(d_msg)


if __name__ == '__main__':
    p = 19
    g = [1,1,1]
    n = 16
    mc = McEliece(p, g ,n)

    e_msg = mc.encrypt(5, mc.publick_key)
    d_msg = mc.decrypt(e_msg, mc.private_key)
    t = 0