# from numpy.lib.function_base import append
from source import *
import numpy as np
from copy import copy
from itertools import zip_longest
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

            
            self.__prime_element = self.__get_prime_element()
            self.__L = self.__gen_L()
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
        def tr(x):
            t = x
            for i in range(self.__m - 1):
                t = mul(t, t, self.__p) ^ x 
            return t

        t = 0
        for i in range(2, 2**self.__m):
            if tr(i) == 1:
                t = i
                break

        return t

    
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
        res = 0
        for i in range(len(g)):
            res ^= self.__pow(x, i + 1)

        return res ^ prime

    
    def __gen_L(self):
        L = []
        for i in range(self.__n):
            L.append(i)

        return L

    
    def __gen_H(self):
        inverse_g = []
        for i in self.__L:
            inverse_g.append(self.__inverse(self.__calc_g(i)))
        
        alpha = [1] * self.__n
        H = []

        for i in range(self.__r):
            H.append([self.__mul(a, inv_g) for a, inv_g in zip(alpha, inverse_g)])
            alpha = [self.__mul(a, l) for a, l in zip(alpha, self.__L)]

        bH = convert_to_binary_matrix(H, self.__m)
        GG = bH

        def find_indx(x, ind):
            for i in range(len(x)):
                if x[i][ind] == 1:
                    return i
            return None

        def swap_columns(x, ind1, ind2):
            for i in range(len(x)):
                x[i][ind1], x[i][ind2] = x[i][ind2], x[i][ind1]
            self.__L[ind1], self.__L[ind2] = self.__L[ind2], self.__L[ind1]

        def find_column_ind(x, ind):
            for i in range(len(x[0])):
                if x[ind][i] == 1:
                    return i

        for i in range(self.__m * self.__t):
            ind = find_indx(GG[i:], i)
            if ind == None:
                ind_c = find_column_ind(GG, i)
                swap_columns(GG, i, ind_c)
                ind = i
            else:
                ind += i

            # ind += i
            max_element = copy(GG[ind])
            # ind_ = self.__k + i
            ind_ = i
            
            GG[ind] = copy(GG[i])
            GG[i] = copy(max_element)

            for j in range(i + 1, self.__m * self.__t):
                if GG[j][ind_] == 1:
                    GG[j] ^= max_element


        mt = self.__m * self.__t

        for i in range(mt - 1):
            for j in range(i  +1, mt):
                if GG[i][j] == 1:
                    GG[i] ^= GG[j]

        
        return GG


    def __gen_G(self):
        return np.concatenate(
            (
                transpose(
                    np.array([i[-self.__k:] for i in self.__H])
                ),
                np.eye(self.__k).astype(np.uint64)
            ), 
            axis=1
        ).astype(np.uint64)


    def __sigma(self, msg):
        def div(divided, divider):
            return div_p(divided, divider, self.__p)

        def mul(multiplier1, multiplier2):
            return mul_p(multiplier1, multiplier2, self.__p)

        def inverse(element):
            ind = 0
            while element[ind] == 0:
                ind += 1
            tmp_e = element[ind:]
            return mod_inverse_p(tmp_e, self.__g, self.__p)
            return mod_inverse_p(element, self.__GC.g, self.__GC.p)


        R = [self.__L[i] for i in range(len(msg)) if msg[i] == 1]
        S_error = [0 for i in range(self.__t)]

        for i in R:
            inv = inverse([1, i])
            # inv = inv[::-1]
            t = 0
            for j in range(len(inv)):
                S_error[j] ^= inv[j]

        # S_error = S_error[::-1]
        if S_error.count(0) == len(S_error):
            return None

        T_error = inverse(S_error)
        
        if T_error == [1, 0]:
            return [1, 0]
        
        if len(T_error) > 1:
            T_error[-2] ^= 1
        else:
            T_error = [1] + T_error

        def get_s():
            e = self.__m * self.__r - 1
            x = [1, 0]
            if e == 0:
                return [1]

            for i in range(e + 1):
                s = copy(x)
                s = div(copy(s), copy(self.__g))[1]
                x = mul(copy(s), copy(x))

            return s

        s = get_s()

        def square_root_p(element):
            element = element[::-1]
            e = 2 ** (self.__m - 1)
            s1 = [self.__pow(element[2 * i], e) for i in range((self.__r - 1) // 2 + 1)][::-1]
            
            s2 = [self.__pow(element[2 * i + 1], e) for i in range(self.__r // 2)][::-1]
            s2 = mul(s2, get_s())

            s = [x1 ^ x2 for x1, x2 in zip_longest(s1[::-1], s2[::-1], fillvalue=0)][::-1]
            return div(list(s), self.__g)[1]

        sr = square_root_p(T_error)

        def mini_euclid(element):
            v = [[0], [1]]
            rem = None
            value = copy(element)
            a = copy(self.__g)

            while (len(value) - 1) > (self.__r // 2) and (len(a) - 1) <= (self.__r // 2):
                t = div(copy(a), copy(value))
                tmp = t[0]
                rem = t[1]

                tmp_mul = mul(copy(v[1]), copy(tmp))
                v[0] = v[0] + [0 for i in range(len(tmp_mul) - len(v[0]))]
                for i in range(len(v[0])):
                    v[0][i] ^= tmp_mul[i]
                v[0], v[1] = v[1], v[0]
                a = copy(value)
                value = copy(rem)

            return [value, v[1]]

        t = mini_euclid(sr)

        r = mul(t[0], t[0])
        v = mul(t[1], t[1])
        v = v + [0]

        sigma = [x1 ^ x2 for x1, x2 in zip_longest(r[::-1], v[::-1], fillvalue=0)][::-1]
        return sigma


    def __errors(self, msg):
        def sigma(_x, pol):
            x = 1
            sum = 0
            pol = pol[::-1]
            for i in pol:
                sum ^= self.__mul(x, i)
                x = self.__mul(x, _x)
            return sum

        sig = self.__sigma(msg)
        if sig == None:
            return np.zeros(self.__n).astype(np.uint64)

        roots = [i for i in self.__L if sigma(i, sig) == 0]
        indx = [self.__L.index(el) for el in roots]
        res = np.zeros(self.__n).astype(np.uint64)
        for i in indx:
            res[i] = 1
        return res


    def error_corrections(self, msg, H = None):
        # знаходимо помилки
        err = self.__errors(msg)

        # прибираємо помилки
        msg = msg ^ err

        return msg


if __name__ == '__main__':
    GoppaCode(19, [1,1,1], 12)