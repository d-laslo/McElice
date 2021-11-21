from numpy.lib.function_base import append
from source import *
import numpy as np
from copy import copy
# from random import randint
import random
import os

class GoppaCode:
    def __init__(self, p, g, n = 0):
        self.__p = p
        self.__n = n

        self.__t = len(g) - 1
        self.__r = self.__t
        self.__m = len(bin(self.__p)[2:]) - 1
        
        if self.__n != 0:
            self.__k = self.__n - self.__m * self.__t
            self.__d = self.__n - self.__k + 1

            self.__L = self.__gen_L()
            self.__prime_element = self.__get_prime_element()
            self.__g = self.__get_g(g)
            self.__H = self.__gen_H()
            self.__G = self.__gen_G()
            # self.__H = self.__tmp_gen_H()
            # t = np.dot(self.__H, np.array(self.__G).T) % 2
            # tt = 0

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
        g.pop(-1)
        g.append(self.__prime_element)
        return g


    def __mul(self, x1, x2):
        return mul(x1, x2, self.__p)


    def __inverse(self, x):
        return mod_inverse(x, self.__p)


    def __pow(self, base, exponenta):
        return pow(base, exponenta, self.__p)


    def __calc_g(self, x):
        g = self.__g[::-1][1:]
        prime = self.__g[-1]
        base = 1
        res = 0
        for i in range(len(g)):
            res ^= self.__pow(x, i + 1) * g[i]
            # res ^= self.__mul(i, base)
            # base = self.__mul(base, x)

        return res ^ prime

    
    def __gen_L(self):
        L = [0, 1]
        for i in range(2, self.__n):
            r = L[-1] << 1
            if len(bin(r)[2:]) > self.__m:
                r ^= self.__p
            L.append(r)

        return L
        # return [4, 8, 3, 6, 12, 11, 5, 10, 7, 14, 15, 13]
        # return [i for i in range(self.__n)]

    
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

        gt_size = len(bH[0])
        gt_size_element = self.__k

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
        # equations = equations
        max_num = 2**gt_size_element - 1 #2**(self.__m * self.__r) - 1
        
        while True:
            # # eye = [] #[2**i for i in range(self.__k)]
            G_T = [0 for i in range(gt_size)]
            # G_T =  [2**i for i in range(self.__k)][::-1] + [0 for i in range(gt_size - self.__k)]
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
                
                #!######################################
                # remove_ind = []
                # if len(eye) > 0 and len(null_indexes) > 1:
                #     for j in range(1, len(null_indexes)):
                #         G_T[null_indexes[j]] = eye.pop(0)
                #         remove_ind.append(j)
                # remove_ind = remove_ind[::-1]
                # for i in remove_ind:
                #     null_indexes.pop(i)     

                #!######################################
                
                for i in null_indexes:
                    random.seed(os.urandom(15))
                    G_T[i] = random.randint(1, max_num)

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
                        tmp = random.randint(1, max_num)

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


            # bed_ind = []
            # # G_T = G_T
            # for ind in range(len(bH)):
            #     sum = 0
            #     for i in bH[ind]:
            #         if i:
            #             sum ^= G_T[i]
            #     if sum !=0:
            #         bed_ind.append(ind)
            # tt = 0


            # GGT = split(G_T, self.__k)
            # tt = np.dot(bH, GGT) %2
            G = binary_transpose(G_T, self.__k)
            
            # tt = split(G, self.__n)
            GG = split(G, self.__n)
            if min(G) == 0:
                continue
            #     if is_transpose:
            #         gg = split(G_T, 12)
            #         return gg
            #     gg = 
            return split(G, self.__n)


    # def __tmp_gen_H(self):
    #     bH = self.__G #convert_to_binary_matrix(self.__H, self.__m)

    #     gt_size = len(bH[0])
    #     gt_size_element = self.__r * self.__m

    #     equations = []
    #     for row in bH:
    #         equations.append(int(''.join([str(em) for em in row]), 2))

    #     # приводимо систему рівнянь до трапецеїдального вигляду
    #     for i in range(len(equations)):
    #         # max_element = max(equations[i:])
    #         ind = index_max_len(equations[i:]) + i
    #         max_element = equations[ind]
    #         len_max_element = len(bin(max_element)) - 2

    #         equations[equations.index(max_element)] = equations[i]
    #         equations[i] = max_element

    #         for j in range(i + 1, len(equations)):
    #             if equations[j] > 0 and len_max_element == (len(bin(equations[j])) - 2):
    #                 equations[j] ^= max_element

    #     tt = len(bin(equations[0])[2:]) - 1
    #     equations = equations[::-1]
    #     # equations = equations
    #     max_num = 2**gt_size_element - 1 #2**(self.__m * self.__r) - 1
        
    #     while True:
    #         # eye = [] #[2**i for i in range(self.__k)]
    #         G_T = [0 for i in range(gt_size)]
    #         # G_T =  [2**i for i in range(self.__k)][::-1] + [0 for i in range(gt_size - self.__k)]
    #         for itr in equations:
    #             indexes = get_indexes(itr)
    #             for i in range(len(indexes)):
    #                 indexes[i] = tt - indexes[i]

    #             # знаходимо індекси всіх незаповнених елементів матриці G
    #             null_indexes = []
    #             for i in indexes:
    #                 if G_T[i] == 0:
    #                     null_indexes.append(i)

    #             if len(null_indexes) == 0:
    #                 continue
                
    #             #!######################################
    #             # remove_ind = []
    #             # if len(eye) > 0 and len(null_indexes) > 1:
    #             #     for j in range(1, len(null_indexes)):
    #             #         G_T[null_indexes[j]] = eye.pop(0)
    #             #         remove_ind.append(j)
    #             # remove_ind = remove_ind[::-1]
    #             # for i in remove_ind:
    #             #     null_indexes.pop(i)     

    #             #!######################################
                
    #             for i in null_indexes:
    #                 random.seed(os.urandom(15))
    #                 G_T[i] = random.randint(1, max_num)

    #             # перший незаповнений елемент є сумою всіх інших елементів у строкі
    #             res = 0
    #             for i in indexes:
    #                 if i == null_indexes[0]:
    #                     continue
    #                 res ^= G_T[i]
    #             G_T[null_indexes[0]] = res if res > 0 else G_T[null_indexes[0]]


    #             # якщо сума всіх елементів, окрім першого незаповненого рівняється 0,
    #             # тоді має існувати, як мінімум два незаповнених елементи, 
    #             # у такому разі достатньо перегенерувати другий незаповнений елемент, 
    #             # і перерахувати перший незаповнений елемент

    #             if G_T[null_indexes[0]] == 0:
    #                 tmp = G_T[null_indexes[1]]
    #                 while tmp == G_T[null_indexes[1]]:
    #                     tmp = random.randint(1, max_num)

    #                 G_T[null_indexes[0]] = G_T[null_indexes[1]] ^ tmp
    #                 G_T[null_indexes[1]] = tmp

    #             # indexes = get_indexes(itr)
    #             # for i in range(len(indexes)):
    #             #     indexes[i] = tt - indexes[i]
    #             sum = 0
    #             for i in indexes:
    #                 sum ^= G_T[i]
    #             ttt = 0

    #         def split(matrix, k):
    #             reslt_matrix = []
    #             _len = k
    #             for row in matrix:
    #                 row = '0' * (_len - len(bin(row)[2:])) + bin(row)[2:]
    #                 reslt_matrix.append([int(i, 2) for i in list(row)])
    #             return reslt_matrix


    #         # bed_ind = []
    #         # # G_T = G_T
    #         # for ind in range(len(bH)):
    #         #     sum = 0
    #         #     for i in bH[ind]:
    #         #         if i:
    #         #             sum ^= G_T[i]
    #         #     if sum !=0:
    #         #         bed_ind.append(ind)
    #         # tt = 0


    #         # GGT = split(G_T, self.__k)
    #         # tt = np.dot(bH, GGT) %2
    #         G = binary_transpose(G_T, gt_size_element)
            
    #         # tt = split(G, self.__n)
    #         GG = split(G, gt_size)
    #         if min(G) == 0:
    #             continue
    #         #     if is_transpose:
    #         #         gg = split(G_T, 12)
    #         #         return gg
    #         #     gg = 
    #         return split(G, gt_size)



if __name__ == '__main__':
   GoppaCode(19, [1,1,1], 12)