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

if __name__ == '__main__':
    pass