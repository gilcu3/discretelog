#! /usr/bin/env python3
from math import sqrt, exp, log, prod

from primefac import isprime
from tqdm import tqdm

from .common import (
    is_Bsmooth,
    smooth_primes,
    phi,
    primitive_root,
    multiplicative_order,
    row_reduce,
)
from .utils import mrange


def producttree(X):
    result = [X]
    while len(X) > 1:
        X = [prod(X[i * 2 : (i + 1) * 2]) for i in range((len(X) + 1) // 2)]
        result.append(X)
    return result


def remaindersusingproducttree(n, T):
    result = [n]
    for t in reversed(T):
        result = [result[i // 2] % t[i] for i in range(len(t))]
    return result


def remainders(n, X):
    return remaindersusingproducttree(n, producttree(X))


def berstein_prec(b):
    ps = smooth_primes(b)
    z = prod(ps)
    return z


# this is probably not worth it in python
def berstein_batch_smooth(b, ns, z):
    pss = []
    zxa = remainders(z, [_[0] for _ in ns])
    for (x, k), zx in zip(ns, zxa):
        e2, y = 2, zx
        while e2 < x and y != 0:
            e2 **= 2
            y = (y * y) % x
        if y == 0:
            _, ps = is_Bsmooth(b, x)
            pss += [(ps, k)]
    return pss


class CongruenceFinder:
    def __init__(self, g, B, p):
        self.g, self.B, self.p = g, B, p
        self.congs = []
        self.k = 0

    def get(self, n, DEBUG=False):
        if DEBUG:
            print(f"Searching {n} congruences")
        batch_size = self.B * 4
        z = berstein_prec(self.B)
        if DEBUG:
            pbar = tqdm(total=n)
        found = 0
        while found < n:
            bk = []
            for _ in range(batch_size):
                self.k += 1
                bk += [(pow(self.g, self.k, self.p), self.k)]
            pss = berstein_batch_smooth(self.B, bk, z)
            found += len(pss)
            self.congs += pss
            if DEBUG:
                pbar.update(len(pss))
        if DEBUG:
            pbar.close()
        bases = list(set(base for c in self.congs for base in c[0]))
        if DEBUG:
            print(f"resulted in bases: {len(bases)}")
        return bases, self.congs


def to_matrices(bases, congruences):
    M = [
        [c[0][base] if base in c[0] else 0 for base in bases] + [c[1]]
        for c in congruences
    ]
    return M


def modMatrix(xx, m):
    x = None
    if isinstance(xx, list):
        x = [None] * (len(xx))
        for i in range(len(xx)):
            if isinstance(xx[i], list):
                x[i] = [None] * len(xx[i])
                for j in range(len(xx[i])):
                    x[i][j] = xx[i][j] % m
            else:
                x[i] = xx[i] % m
    return x


def evaluate(eq, dlogs, q):
    return sum([dlogs[term] * exp for term, exp in eq.items()]) % q


def check_congruences(congruences, dlogs, q):
    passed = True
    for c in congruences:
        if evaluate(c[0], dlogs, q) != c[1] % q:
            passed = False
            break
    return passed


def check_dlogs(g, p, q, exponents, bases):
    passed = True
    o = phi(p)
    for exponent, base in zip(exponents, bases):
        if pow(g, exponent * o // q, p) != pow(base, o // q, p):
            passed = False
            break
    return passed


def msolve_prime(M, q, DEBUG=False):
    if DEBUG:
        print(f"solving linear system {len(M)}x{len(M[0])}")
    n = len(M[0]) - 1
    m = modMatrix(M, q)
    if not row_reduce(m, q):
        return None
    return [m[i][n] for i in range(n)]


def dlog_prime(b, h, p, DEBUG=False):
    if p <= 10**5:
        raise AssertionError(f"{p} is too small")
    q = multiplicative_order(b, p)
    assert isprime(q)
    o = phi(p)
    g = primitive_root(p)
    B = int(4 * exp(sqrt(log(p) * log(log(p))) / 2)) + 10
    if DEBUG:
        print("p: {}, b: {}, h: {}, B: {}".format(p, b, h, B))
    congs = []
    cf = CongruenceFinder(g, B, p)
    while True:
        crels = len(smooth_primes(B)) + 1 if len(congs) == 0 else 50
        bases, congs = cf.get(crels, DEBUG)
        m = to_matrices(bases, congs)
        exponents = msolve_prime(m, q)
        if exponents is not None:
            break
    if DEBUG:
        print("done.")

    dlogs = {base: exp for (base, exp) in zip(bases, exponents)}

    assert check_congruences(congs, dlogs, q)
    assert check_dlogs(g, p, q, exponents, bases)

    def find_log_g(x):
        if DEBUG:
            print(f"searching for k such that {x}*g^-k is B-smooth.")
        for k in mrange(1, 10**7, 1, DEBUG):
            c = is_Bsmooth(B, x * pow(g, -k, p) % p)
            if c[0] and set(c[1].keys()).issubset(set(dlogs.keys())):
                if DEBUG:
                    print("found k = {}".format(k))
                break

        xe = (evaluate(c[1], dlogs, q) + k) % q
        if DEBUG:
            print(f"Found partial dlog: {g}^{xe}={x} mod {p} mod {q}")
        assert pow(g, xe * o // q, p) == pow(x, o // q, p)
        return xe

    be = find_log_g(b)
    he = find_log_g(h)

    sol = he * pow(be, -1, q) % q
    if pow(b, sol, p) == h:
        if DEBUG:
            print("{}^{} = {} (mod {}) holds!".format(b, sol, h, p))
            print("DLP solution: {}".format(sol))
        return sol
