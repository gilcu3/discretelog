#! /usr/bin/env python
from .utils import mrange


def to_sp(mat):
    n, m = len(mat), len(mat[0])
    nmat = [[] for _ in range(n)]
    for i in range(n):
        for j, v in enumerate(mat[i]):
            if v != 0:
                nmat[i].append((j, v))
    return nmat, m


def from_sp(mat, m):
    n = len(mat)
    nmat = [[0] * m for _ in range(n)]
    for i in range(n):
        for j, v in mat[i]:
            nmat[i][j] = v
    return nmat


def sptr(mat, m):
    nmat = [[] for _ in range(m)]
    n = len(mat)
    for i in range(n):
        for j, v in mat[i]:
            nmat[j].append((i, v))
    return nmat


def spmul(mat, vec, q):
    n = len(mat)
    nvec = [0] * n
    for i in range(n):
        nvec[i] = sum(v * vec[j] for j, v in mat[i]) % q
    return nvec


def dot(vec1, vec2):
    return sum(v1 * v2 for v1, v2 in zip(vec1, vec2))


def smul(c, vec, q):
    return [c * v for v in vec]


def vadd(vec1, vec2):
    return [v1 + v2 for v1, v2 in zip(vec1, vec2)]


def vsub(vec1, vec2):
    return [v1 - v2 for v1, v2 in zip(vec1, vec2)]


def vmod(vec, q):
    return [c % q for c in vec]


def mod(a, q):
    if abs(a) >= q // 2:
        a = a % q
        if a > q // 2:
            a -= q
    return a


def block_lanczos(mat, d, b, q, DEBUG=False):
    m, mt = mat, sptr(mat, d)
    v0 = spmul(mt, b, q)

    def mulA(x):
        return spmul(mt, spmul(m, x, q), q)
    x = [0] * d
    denp = None
    v2, v1 = None, v0
    for i in mrange(1, d + 1, DEBUG=DEBUG):
        Avi = mulA(v1)
        num = dot(Avi, Avi) % q
        den1 = dot(v1, Avi) % q
        if den1 % q == 0:
            # it must finish exactly on the d-th iteration
            return None
        den1 = pow(den1, -1, q)
        ci = smul(num * den1 % q, v1, q)
        numx = dot(v1, v0) % q
        x = vmod(vadd(x, smul(numx * den1 % q, v1, q)), q)
        vi = vsub(Avi, ci)
        if i > 1:
            num = dot(v1, Avi) % q
            ci1 = smul(num * denp % q, v2, q)
            vi = vsub(vi, ci1)
        denp = den1
        v2, v1 = v1, vmod(vi, q)
    if not all(c % q == 0 for c in v1):
        return None
    # assert b == vmod(spmul(mat, x), q)
    return vmod(x, q)


if __name__ == '__main__':
    m, d = to_sp([[1, 2], [3, 2], [0, 1]])
    q = 2 ** (2 ** 3) + 1
    b = [-1, 1, -1]
    x = block_lanczos(m, d, b, q)
    print(x)
