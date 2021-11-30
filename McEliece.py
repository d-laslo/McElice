from GoppaCode import GoppaCode as GC
from source import *
from random import sample
import numpy as np

class McEliece:
    def __init__(self, p, g, n = 0):
        self.__GC_flag = False
        if n != 0:
            self.__GC_flag = True
        self.__GC = GC(p, g, n)
        self.__G = []
        self.__H = []

        self.__S = []
        self.__P = []
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
            t = np.linalg.det(S) % 2

        S = np.array([
            [1, 1, 0, 0],
            [0, 0, 1, 1],
            [0, 0, 0, 1],
            [0, 1, 1, 0]
        ], dtype=np.uint64)
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
        return [self.__S, self.__G, self.__P]


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


    def encrypt(self, msg, public_key):
        Gh, t = public_key
        k = len(Gh)
        n = len(Gh[0])

        err = np.zeros(n, np.uint64)
        for i in sample(range(n), t):
            err[i] = 1

        bmsg = self.__convert_to_binary_vector(msg, k)
        return self.__convert_to_num(self.__dot(bmsg, Gh) ^ err)


    def decrypt(self, msg, private_key):
        S, G, P = private_key
        
        # домножуємо на оберену матрицю P
        bmsg = self.__convert_to_binary_vector(msg, self.__GC.n)
        bmsg = self.__dot(bmsg, self.__calc_inv_P(P))

            
        bmsg = self.__GC.error_corrections(bmsg)
        bmsg = bmsg[-self.__GC.k:]

        # домножуємо на оберену матрицю P
        d_msg = self.__dot(bmsg, self.__calc_inv_S(S))

        return self.__convert_to_num(d_msg)


if __name__ == '__main__':
    p = 19
    g = [1,1,1]
    n = 12
    mc = McEliece(p, g ,n)

    e_msg = mc.encrypt(2, mc.publick_key)
    d_msg = mc.decrypt(e_msg, mc.private_key)
    print(d_msg)
    t = 0