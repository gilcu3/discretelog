#! /usr/bin/env python3 -u

from math import isqrt, pi

from primefac import isprime

from .common import order, factor, solve_congruence, crt
from .utils import mrange, mexec
from .linear_sieve_index_calculus import linear_sieve_dlog


def baby_steps_giant_steps(g, v, q):
    po = order(g, q)
    n = isqrt(po) + 1
    an = pow(g, n, q)
    vals = {}
    cur = v
    for w in range(n + 1):
        if cur in vals:
            break
        vals[cur] = w
        cur = cur * g % q

    cur = 1
    for p in range(1, n + 1):
        cur = cur * an % q
        if cur in vals:
            d = n * p - vals[cur]
            return d % po
    return None


def pollard_rho_dlg(g, v, q, DEBUG=False):
    if DEBUG:
        print(f'pollard_rho_dlg: g={g} v={v} q={q}')
    if g == v:
        return 1
    n = order(g, q)
    new_xab = [
        lambda x, a, b: (x * x % q, a * 2 % n, b * 2 % n),
        lambda x, a, b: (x * g % q, (a + 1) % n, b),
        lambda x, a, b: (x * v % q, a, (b + 1) % n)
    ]
    sn = 4 * isqrt(int(pi * n / 8))
    x, a, b = 1, 0, 0
    X, A, B = x, a, b
    for _ in mrange(1, sn, DEBUG=DEBUG):
        x, a, b = new_xab[x % 3](x, a, b)
        X, A, B = new_xab[X % 3](X, A, B)
        X, A, B = new_xab[X % 3](X, A, B)
        if x == X:
            break

    for d in solve_congruence((B - b) % n, (a - A) % n, n):
        if pow(g, d, q) == v:
            if DEBUG:
                print(f'discrete log: d={d}')
            return d

    if DEBUG:
        # Consequence of this bad case
        # 345259615358946849545 720882080607834014753 1097036420717452490927
        print('Falling back to baby steps giant steps')
    d = baby_steps_giant_steps(g, v, q)
    assert pow(g, d, q) == v, f'pollard_rho_dlg failed g={g} v={v} q={q}'
    return d


def dlog_cado_nfs(p, g, v, q, DEBUG=False):
    ds = mexec('cado-nfs -dlp -ell %d target=%d,%d %d' % (p, g, v, q), DEBUG)
    ds = ds.split('\n')[-1]
    dg, dv = tuple(map(int, ds.split(',')))
    d = dv * pow(dg, -1, p) % p
    return d


def dlog_prime_order(g, v, q, p, DEBUG=False):
    if p <= 10 ** 10:
        d = baby_steps_giant_steps(g, v, q)
    elif p <= 10 ** 14:
        # should work for powers of a prime as well
        d = pollard_rho_dlg(g, v, q, DEBUG)
    elif q <= 10 ** 35:
        d = linear_sieve_dlog(q, g, v, p, DEBUG=DEBUG)
    else:
        assert isprime(q), f'modulus q={q} is not prime, cado-nfs requires it'

        # this should never fail
        assert isprime(p), f'order p={p} is not prime, cado-nfs requires it'

        d = dlog_cado_nfs(p, g, v, q, DEBUG)

    assert pow(g, d, q) == v % q, 'discrete log cyclic %d is incorrect' % d

    return d


def dlog_prime_power_order(g, h, q, p, e, DEBUG=False):
    if e == 1:
        return dlog_prime_order(g, h, q, p, DEBUG)
    else:
        pe = p ** e
        x = 0
        pk = pe // p
        gamma = pow(g, pk, q)
        for _ in range(e):
            hk = pow(pow(g, -x, q) * h, pk, q)
            dk = dlog_prime_order(gamma, hk, q, p, DEBUG)
            x = (x + pe // pk // p * dk) % pe
            pk //= p
        assert pow(g, x, q) == h % q, f'dlog prime power order x={x} failed'
        return x


def pohlig_hellman(g, v, q, DEBUG=False):

    d = order(g, q)
    fd = factor(d, DEBUG)
    if len(fd) == 1:
        p, e = list(fd.items())[0]
        x = dlog_prime_power_order(g, v, q, p, e, DEBUG)
    else:
        congs = []
        for p, e in fd.items():
            pi = p ** e
            gi = pow(g, d // pi, q)
            vi = pow(v, d // pi, q)
            xi = dlog_prime_power_order(gi, vi, q, p, e, DEBUG)
            congs += [(xi, pi)]
        x, _ = crt(congs)
    assert pow(g, x, q) == v, \
        f'pohlig hellman failed: g={g} v={v} d={d} x={x} q={q}'
    return x


def dlog_prime_power(g, v, p, e, DEBUG=False):

    q = p ** e
    o = order(g, q)
    assert pow(v, o, q) == 1, \
        '%d is not in the subgroup of %d mod %d' % (v, g, q)

    if q <= 10 ** 10:

        d = baby_steps_giant_steps(g, v, q)

    else:

        d = pohlig_hellman(g, v, q, DEBUG)

    assert pow(g, d, q) == v, \
        'discrete log prime power %d is incorrect' % d
    return d, o


def dlog_non_relative_prime(g, v, p, e):
    q = p ** e
    if g == 0 and v == 0:
        return 1, 1

    gd = 1
    for d in range(e + 1):
        if gd == v:
            if d == e or gd == 0:
                return d, 1
            else:
                return d, 0
        gd = gd * g % q
    return None, None


def dlog(g, v, n, DEBUG=False):
    fn = factor(n, DEBUG)
    congs = []
    c, ct = None, None

    for p, e in fn.items():
        if g % p == 0:
            q = p ** e
            gi, vi = g % q, v % q
            xi, ti = dlog_non_relative_prime(gi, vi, p, e)
            if ct is None:
                c, ct = xi, ti
            else:
                if ti == 0:
                    if ct == 0 and c != xi:
                        return None
                    elif ct == 1 and c > xi:
                        return None
                    else:
                        c, ct = xi, ti
                elif ti == 1:
                    if ct == 0:
                        if xi > c:
                            return None
                    elif ct == 1:
                        c = max(c, xi)
    if ct == 0:
        if pow(g, c, n) == v:
            return c, 0
        else:
            return None

    for p, e in fn.items():
        if g % p != 0:
            pi = p ** e
            gi = g % pi
            vi = v % pi
            if gi == 1:
                if vi != 1:
                    return None
            else:
                xi, ti = dlog_prime_power(gi, vi, p, e, DEBUG)
                congs += [(xi, ti)]

    d, m = crt(congs)
    dt = m

    if ct is not None:
        assert ct == 1
        if d < c:
            d = d + m * ((c - d + m - 1) // m)

    assert pow(g, d, n) == v, 'discrete log %d is incorrect' % d
    return d, dt
