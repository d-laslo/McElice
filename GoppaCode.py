from numpy.lib.function_base import append
from source import *
import numpy as np
from copy import copy
from random import randint

class GoppaCode:
    def __init__(self, p, g, n = 0):
        self.__p = p
        self.__n = n

        self.__t = (len(g) - 1) // 2
        self.__r = len(g) - 1
        self.__m = len(bin(self.__p)[2:]) - 1
        
        if self.__n != 0:
            self.__k = self.__n - self.__r
            self.__d = self.__n - self.__k + 1

            self.__L = self.__gen_L()
            self.__prime_element = self.__get_prime_element()
            self.__g = self.__get_g(g)
            self.__H = self.__gen_H()
            self.__G = self.__gen_G()

    @property
    def H(self):
        return self.__H


    @property
    def G(self):
        return self.__G
    

    @property
    def n(self):
        return self.__n


    @property
    def p(self):
        return self.__p


    @property
    def t(self):
        return self.__t


    @property
    def r(self):
        return self.__r


    @property
    def m(self):
        return self.__m


    @property
    def k(self):
        return self.__k


    @property
    def d(self):
        return self.__d

    
    @property
    def L(self):
        return self.__L


    @property
    def g(self):
        return self.__g

    
    def __get_prime_element(self):
        return 8

    
    def __get_g(self, g):
        g[-1] = self.__prime_element
        return g


    def __mul(self, x1, x2):
        return mul(x1, x2, self.__p)


    def __inverse(self, x):
        return mod_inverse(x, self.__p)


    def __pow(self, base, exponenta):
        return pow(base, exponenta, self.__p)


    def __calc_g(self, x):
        g = self.__g[::-1][1:]
        base = 1
        res = 0
        for i in range(len(g)):
            res ^= self.__pow(x, i + 1) * g[i]
            # res ^= self.__mul(i, base)
            # base = self.__mul(base, x)

        return res ^ self.__prime_element

    
    def __gen_L(self):
        return [4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13]
        return [i for i in range(self.__n)]

    
    def __gen_H(self):
        inverse_g = []
        for i in self.__L:
            inverse_g.append(self.__inverse(self.__calc_g(i)))
        
        alpha = [1] * self.__n
        H = []

        for i in range(self.__r):
            H.append([self.__mul(a, inv_g) for a, inv_g in zip(alpha, inverse_g)])
            alpha = [self.__mul(a, l) for a, l in zip(alpha, self.__L)]

        return H


    def __gen_G(self):
        bH = convert_to_binary_matrix(self.__H, self.__m)

        # is_transpose = False
        

        # if len(bH) > len(bH[0]):
        #     bH = transpose(bH)
        #     is_transpose = True
        #     # gt_size = self.__k
        #     # gt_size_element = self.__n

        gt_size = len(bH[0])
        gt_size_element = self.__k

        # for row in H_T:
        #     brow = ['0' * (self.__m - len(bin(el)[2:])) + bin(el)[2:] for el in row]
        #     bH_T.append(list(''.join(brow)))

        # bH = list(np.array(bH_T).transpose())

        equations = []
        for row in bH:
            equations.append(int(''.join([str(em) for em in row]), 2))

        # приводимо систему рівнянь до трапецеїдального вигляду
        for i in range(len(equations)):
            # max_element = max(equations[i:])
            ind = index_max_len(equations[i:]) + i
            max_element = equations[ind]
            len_max_element = len(bin(max_element)) - 2

            equations[equations.index(max_element)] = equations[i]
            equations[i] = max_element

            for j in range(i + 1, len(equations)):
                if equations[j] > 0 and len_max_element == (len(bin(equations[j])) - 2):
                    equations[j] ^= max_element

        tt = len(bin(equations[0])[2:]) - 1
        equations = equations[::-1]
        max_num = 2**gt_size_element - 1 #2**(self.__m * self.__r) - 1
        while True:
            G_T = [0 for i in range(gt_size)] 
            for itr in equations:
                indexes = get_indexes(itr)
                for i in range(len(indexes)):
                    indexes[i] = tt - indexes[i]

                # знаходимо індекси всіх незаповнених елементів матриці G
                null_indexes = []
                for i in indexes:
                    if G_T[i] == 0:
                        null_indexes.append(i)

                if len(null_indexes) == 0:
                    continue

                # if len(null_indexes) == 1:
                #     G_T[null_indexes[0]] = -1
                #     continue
                
                for i in null_indexes:
                    G_T[i] = randint(1, max_num)

                # перший незаповнений елемент є сумою всіх інших елементів у строкі
                res = 0
                for i in indexes:
                    if i == null_indexes[0]:
                        continue
                    res ^= G_T[i]
                G_T[null_indexes[0]] = res if res > 0 else G_T[null_indexes[0]]


                # якщо сума всіх елементів, окрім першого незаповненого рівняється 0,
                # тоді має існувати, як мінімум два незаповнених елементи, 
                # у такому разі достатньо перегенерувати другий незаповнений елемент, 
                # і перерахувати перший незаповнений елемент

                if G_T[null_indexes[0]] == 0:
                    tmp = G_T[null_indexes[1]]
                    while tmp == G_T[null_indexes[1]]:
                        tmp = randint(1, max_num)

                    G_T[null_indexes[0]] = G_T[null_indexes[1]] ^ tmp
                    G_T[null_indexes[1]] = tmp

                # indexes = get_indexes(itr)
                # for i in range(len(indexes)):
                #     indexes[i] = tt - indexes[i]
                # sum = 0
                # for i in indexes:
                #     sum ^= G_T[i]
                # ttt = 0

            def split(matrix, k):
                reslt_matrix = []
                _len = k
                for row in matrix:
                    row = '0' * (_len - len(bin(row)[2:])) + bin(row)[2:]
                    reslt_matrix.append([int(i, 2) for i in list(row)])
                return reslt_matrix

            G = binary_transpose(G_T, self.__k)
            # if min(G) > 0:
            #     if is_transpose:
            #         gg = split(G_T, 12)
            #         return gg
            #     gg = 
            return split(G, self.__n)



if __name__ == '__main__':
   GoppaCode(12, 19, [1, 0, 1, 0, 1])
#    GoppaCode(12, 19, [1, 1, 1])