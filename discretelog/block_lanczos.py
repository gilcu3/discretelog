#! /usr/bin/env python

from functools import partial
import sys

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


def spmax(mat):
    mx = 0
    for i in range(len(mat)):
        c1, c2 = 0, 0
        for _, v in mat[i]:
            if v > 0:
                c1 += v
            else:
                c2 -= v
        mx = max(mx, max(c1, c2))
    return mx


def spmul2(mx, mat, vec, q):
    n = len(mat)
    nvec = [0] * n
    xvec = vec[:]
    bt = 63 - mx.bit_length() - 1
    mask = (1 << bt) - 1
    for b in range(0, q.bit_length() + 1, bt):
        cvec = [v & mask for v in xvec]
        for i in range(n):
            a = 0
            for j, v in mat[i]:
                a += cvec[j] * v
                # assert abs(a) <= 1 << 63
            nvec[i] += a << b
        for i in range(len(xvec)):
            xvec[i] >>= bt
    return vmod(nvec, q)


def spmul(mat, vec, q):
    n = len(mat)
    nvec = [0] * n
    for i in range(n):
        a = 0
        for j, v in mat[i]:
            a += vec[j] * v
        nvec[i] = a
        # for some reason this is slower
        # nvec[i] = sum([v * vec[j] for j, v in mat[i]])
    return vmod(nvec, q)


def dot(vec1, vec2, q):
    return sum([v1 * v2 for v1, v2 in zip(vec1, vec2)]) % q


def smul(c, vec, q):
    return [c * v % q for v in vec]


def vadd(vec1, vec2):
    return [v1 + v2 for v1, v2 in zip(vec1, vec2)]


def vsub(vec1, vec2):
    return [v1 - v2 for v1, v2 in zip(vec1, vec2)]


def vmod(vec, q):
    return [c % q for c in vec]


def matmul_choice(q, m, DEBUG=False):
    # This function chooses the best matrix multiplication routine
    # In my tests spmul2 with integer unpacked works better using pypy3 for
    # bigger numbers, which probably should not be the case

    is_pypy = "PyPy" in sys.version
    if is_pypy:
        mx = spmax(m)
        if q.bit_length() + mx.bit_length() > 70 and mx.bit_length() <= 15:
            if DEBUG:
                print(f"Using packed matmul for m with mx={mx.bit_length()}")
            return partial(spmul2, mx)
    return spmul


def block_lanczos(mat, d, b, q, DEBUG=False):
    m, mt = mat, sptr(mat, d)

    # This is a hack to achieve better matrix multiplication performance
    rmul = matmul_choice(q, m, DEBUG)
    rmult = matmul_choice(q, mt, DEBUG)

    v0 = rmult(mt, b, q)

    def mulA(x):
        return rmult(mt, rmul(m, x, q), q)

    x = [0] * d
    denp = None
    v2, v1 = None, v0
    for i in mrange(1, d + 1, DEBUG=DEBUG):
        Avi = mulA(v1)
        num = dot(Avi, Avi, q)
        den1 = dot(v1, Avi, q)
        if den1 % q == 0:
            # it must finish exactly on the d-th iteration
            return None
        den1i = pow(den1, -1, q)
        ci = smul(num * den1i % q, v1, q)
        numx = dot(v1, v0, q)
        x = vmod(vadd(x, smul(numx * den1i % q, v1, q)), q)
        vi = vsub(Avi, ci)
        if i > 1:
            num = den1
            ci1 = smul(num * denp % q, v2, q)
            vi = vsub(vi, ci1)
        denp = den1i
        v2, v1 = v1, vmod(vi, q)
    if not all(c % q == 0 for c in v1):
        return None
    # assert b == vmod(spmul(mat, x), q)
    return vmod(x, q)
