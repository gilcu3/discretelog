#! /usr/bin/env python3

from functools import cache
from math import gcd, log2
import random

from primefac import primefac, isprime

from .utils import mexec, mrange


def perfect_power(n):
    for e in range(2, int(log2(n)) + 1):
        r = round(n ** (1 / e))
        if r ** e == n:
            return r, e
    return n, 1


def solve_congruence(a, b, n):
    g = gcd(a, n)
    if b % g != 0:
        return
    ap, bp, np = a // g, b // g, n // g
    for s in range(bp * pow(ap, -1, np) % np, n, n // g):
        yield s


def egcd(a, b):
    if b == 0:
        return (1, 0, a)
    q = egcd(b, a % b)
    return (q[1], q[0] - a // b * q[1], q[2])


def crt(rm):

    def combine(c1, c2):
        a1, m1 = c1
        a2, m2 = c2
        a, b, c = m1, a2 - a1, m2
        g = gcd(a, b, c)
        a, b, c = [i//g for i in [a, b, c]]
        if a != 1:
            inv_a, _, g = egcd(a, c)
            if g != 1:
                return None
            b *= inv_a
        a, m = a1 + m1*b, m1*c
        return a, m

    rv = (0, 1)
    n, m = 0, 1
    for rmi in rm:
        rv = combine(rv, rmi)
        if rv is None:
            break
        n, m = rv
        n = n % m
    else:
        return n, m


def factor_cado_nfs(n, DEBUG=False):
    ds = mexec('cado-nfs %d' % n, DEBUG)
    ds = ds.split('\n')[-1]
    ps = list(map(int, ds.split()))
    nn = n
    pss = []
    for p in ps:
        while nn % p == 0:
            nn //= p
            pss += [p]
    return pss


def factor_yafu(n, pretest=False, DEBUG=False):
    pretest_opt = ' -pretest 20' if pretest else ''
    f = mexec('yafu%s %d' % (pretest_opt, n), DEBUG=DEBUG)
    op = False
    ps = []
    for line in f.split('\n'):
        if line.startswith('***factors found***'):
            op = True
        if op and line.startswith('P'):
            p = int(line.split()[-1])
            ps += [p]
        if op and line.startswith('ans'):
            break
    return ps


@cache
def factor(n, DEBUG=False):
    if DEBUG:
        print(f'factor: n={n}')
    if isprime(n):
        f = {n: 1}
    else:
        f = {}
        if n <= 10 ** 20:
            ps = list(primefac(n, verbose=DEBUG))
        else:
            ps = []
            try:
                for p in primefac(n, verbose=True, trial=10**4, rho=2 * 10**5,
                                  methods=tuple()):
                    if isprime(p):
                        n //= p
                        ps += [p]
                    else:
                        break
            except AssertionError:
                pass
            if n > 1:
                for p in factor_yafu(n, pretest=True, DEBUG=DEBUG):
                    if isprime(p):
                        n //= p
                        ps += [p]
                    else:
                        break
                if n > 1:
                    # This is a hard to factor number
                    if n <= 10 ** 100:
                        ps += factor_yafu(n, DEBUG=DEBUG)
                    else:
                        ps += factor_cado_nfs(n, DEBUG=DEBUG)
        for p in ps:
            if p in f:
                f[p] += 1
            else:
                f[p] = 1
    if DEBUG:
        print(f'f={f}')
    return f


@cache
def phi(q):
    ans = 1
    for p, e in factor(q).items():
        ans *= p ** (e - 1) * (p - 1)
    return ans


@cache
def order(g, q):
    assert gcd(g, q) == 1
    d = phi(q)
    qs = factor(d).keys()
    while True:
        found = True
        assert pow(g, d, q) == 1
        for p in qs:
            if d % p == 0 and pow(g, d // p, q) == 1:
                d //= p
                found = False
                break
        if found:
            break
    assert pow(g, d, q) == 1, f'order failed: g={g} d={d} q={q}'
    return d


@cache
def multiplicative_order(g, q):
    assert gcd(g, q) == 1
    d = phi(q)
    qs = factor(d).keys()
    while True:
        found = True
        assert pow(g, d, q) == 1
        for p in qs:
            if d % p == 0 and pow(g, d // p, q) == 1:
                d //= p
                found = False
                break
        if found:
            break
    assert pow(g, d, q) == 1, f'order failed: g={g} d={d} q={q}'
    return d


@cache
def primitive_root(q):
    o = phi(q)
    for g in range(2, q):
        if multiplicative_order(g, q) == o:
            return g


@cache
def smooth_primes(b):
    primes = [True] * b
    pp = [2]
    for i in range(2, b, 2):
        primes[i] = False
    for i in range(3, b, 2):
        if primes[i]:
            pp += [i]
            for j in range(i ** 2, b, 2 * i):
                primes[j] = False
    return pp


def is_Bsmooth(b, n):
    ps = {}
    for p in smooth_primes(b):
        if n > 1 and p**2 > n:
            if n < b:
                ps[n] = 1
                return True, ps
            else:
                return False, None
        e = 0
        while n % p == 0:
            e += 1
            n //= p
        if e > 0:
            ps[p] = e
    return n == 1, ps


def random_sophie_germain_prime(d, frandom=random):
    while True:
        p = frandom.randint(10 ** d, 10 ** (d + 1))
        if isprime(p) and isprime(2 * p + 1):
            return 2 * p + 1


def random_prime(d, frandom=random):
    while True:
        p = frandom.randint(10 ** d, 10 ** (d + 1))
        if isprime(p):
            return p


def row_reduce(M, q, DEBUG=False):
    q0, _ = perfect_power(q)
    assert isprime(q0)
    n, m = len(M), len(M[0])
    cj = -1
    perfect = True
    for i in mrange(0, m - 1, 1, DEBUG):
        cj += 1
        pivot = None
        for j in range(i, n):
            if M[j][cj] % q0 != 0:
                pivot = j
                break
        if pivot is None:
            perfect = False
            continue
        M[i], M[pivot] = M[pivot], M[i]
        iv = pow(M[i][cj], -1, q)
        if iv != 1:
            for k in range(cj, m):
                M[i][k] = M[i][k] * iv % q
        for j in range(n):
            if j != i:
                if M[j][cj] != 0:
                    for k in range(cj + 1, m):
                        M[j][k] = (M[j][k] - M[i][k] * M[j][cj]) % q
                    M[j][cj] = 0
    return perfect
