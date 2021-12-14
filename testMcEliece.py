from McEliece import McEliece
from source import *
from random import randint
import cProfile

def test_vector1():
    p = convert2num('z^12 + z^3 + 1')
    g = convert2vector('y^64 + y^3 + y + 1')
    n = 3488
    mc = McEliece(p, g ,n)

    e_msg = mc.encrypt(5, mc.publick_key)
    d_msg = mc.decrypt(e_msg, mc.private_key)

def custom_test(times = 1000):
    p = 19
    g = [1,1,1]
    n = 12
    mc = McEliece(p, g ,n)

    res = []

    for i in range(times):
        rr = randint(1, 15)
        e_msg = mc.encrypt(rr, mc.publick_key)
        d_msg = mc.decrypt(e_msg, mc.private_key)
        res.append(rr == d_msg)
    print(res.count(True))

if __name__ == '__main__':
    custom_test(times = 10000);
    # test_vector1()