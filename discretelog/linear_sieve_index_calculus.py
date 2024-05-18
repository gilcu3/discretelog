#! /usr/bin/env python

from math import isqrt, log, exp, ceil
from os import cpu_count

from primefac import isprime

from .common import primitive_root, smooth_primes, factor, is_Bsmooth, \
    row_reduce, order
from .utils import mrange, parallel_for_balanced
from .block_lanczos import to_sp, block_lanczos


def L(p, alpha, c):
    return int(exp(c * log(p)**alpha * log(log(p)) ** (1 - alpha)))


def inverse_seq(x, q):
    qpow = 1
    y, yq = 0, pow(x % q, -1, q)
    while True:
        s = (x % (qpow * q) * y - 1) // (qpow)
        k = (-s * yq) % q
        y += qpow * k
        qpow *= q
        yield y, qpow


def well_connected(k, graph, degc):
    que = []
    seen = [False] * k
    for v in range(k):
        if degc[v] == 1:
            que += [v]
            seen[v] = True
    front = 0
    while front < len(que):
        cv = que[front]
        front += 1
        for v in graph[cv]:
            if not seen[v]:
                degc[v] -= 1
                if degc[v] == 1:
                    seen[v] = True
                    que += [v]


def structured_gaussian_elimination(rels, relsex, k):
    cols = [set() for _ in range(k)]
    kp = len(rels[0])
    n = len(rels)
    colsw = [set() for _ in range(n)]
    singlerow = []
    doublerow = []
    for i, crex in enumerate(relsex):
        if len(crex) == 1:
            singlerow += [i]
        if len(crex) == 2:
            doublerow += [i]
        for c in crex:
            cols[c].add(i)
    for c in range(k):
        w = len(cols[c])
        colsw[w].add(c)

    assert len(colsw[0]) == 0

    def reduce_row(c, orig, new):
        relsex[new].remove(c)
        if len(relsex[new]) == 1:
            singlerow.append(new)
        for i in range(kp):
            rels[new][i] -= rels[orig][i]

    marked = [False] * n
    rn = n
    changed2 = True
    extrarels = 1.1
    while changed2:
        changed = True
        while changed:
            changed = False
            # delete single rows
            tmp = singlerow[:]
            if len(singlerow) > 0:
                singlerow.clear()
                changed = True
            for r in tmp:
                if len(relsex[r]) == 1 and not marked[r]:
                    c = relsex[r][0]
                    if len(cols[c]) > 1:
                        for i in cols[c]:
                            if i != r and not marked[i]:
                                reduce_row(c, r, i)
                        colsw[len(cols[c])].remove(c)
                        cols[c] = set([r])
                        colsw[1].add(c)
            # delete single columns
            torm = []
            for cr in colsw[1]:
                assert len(cols[cr]) == 1
                for r in cols[cr]:
                    if not marked[r]:
                        torm += [r]
                        marked[r] = True
                        rn -= 1
            if len(colsw[1]) > 0:
                changed = True
            for r in torm:
                for c in relsex[r]:
                    cols[c].remove(r)
                    colsw[len(cols[c]) + 1].remove(c)
                    colsw[len(cols[c])].add(c)
        changed2 = False
        for c in range(k):
            assert c in colsw[len(cols[c])], (c, cols[c], colsw[len(cols[c])])
        while rn > extrarels + k - len(colsw[0]) + kp and len(doublerow) > 0:
            r = doublerow.pop()
            if len(relsex[r]) == 2 and not marked[r]:
                changed2 = True
                marked[r] = True
                rn -= 1
                for c in relsex[r]:
                    cols[c].remove(r)
                    colsw[len(cols[c]) + 1].remove(c)
                    colsw[len(cols[c])].add(c)
    # this case handles the case that too many relations remain
    if rn > extrarels * (k - len(colsw[0]) + kp):
        assert len(doublerow) == 0
        assert len(colsw[0]) == k
        rw = []
        for i in range(n):
            if not marked[i]:
                w = len([v for v in rels[i] if v != 0])
                rw += [(w, i)]
        rw.sort()
        for _w, i in reversed(rw):
            rn -= 1
            marked[i] = True
            if rn <= kp * extrarels:
                break
    return lambda i: not marked[i]


def mod_matrix(m, p):
    for i in range(len(m)):
        for j in range(len(m[i])):
            m[i][j] %= p
            if m[i][j] > p // 2:
                m[i][j] -= p


def rels_filter(rels, relsex, k):
    degc = [0] * k
    graph = [[] for _ in range(k)]
    for crex in relsex:
        for c in crex:
            degc[c] += 1
        if len(crex) == 2:
            u, v = crex
            graph[u] += [v]
            graph[v] += [u]
    well_connected(k, graph, degc)
    n = len(rels)
    marked = [False] * n
    for i, (_cr, crex) in enumerate(zip(rels, relsex)):
        pos = False
        for c in crex:
            if degc[c] <= 1:
                pos = True
        marked[i] = pos
    return lambda i: not marked[i]


def rels_matrix(rels, relsex, k, rf, DEBUG=False):
    m = 0
    kmap = [None] * k
    nrels = []
    for i, crex in enumerate(relsex):
        if rf(i):
            for c in crex:
                if kmap[c] is None:
                    kmap[c] = m
                    m += 1
    if DEBUG:
        sparsity = 0
        mx = 0
    for i, (cr, crex) in enumerate(zip(rels, relsex)):
        if rf(i):
            cur = [0] * m + cr
            if DEBUG:
                sparsity += len(crex) + len([c for c in cr if c != 0])
                mx = max(mx, max([abs(c) for c in cr]))
            if len(crex) == 2:
                u, v = crex
                cur[kmap[u]] = -2
                cur[kmap[v]] = -2
            elif len(crex) == 1:
                v = crex[0]
                cur[kmap[v]] = -2
            nrels += [cur]
    if DEBUG:
        print(f'system sparsity: {sparsity / len(nrels)}')
        print(f'max absolute value: {mx}')
    return nrels, kmap


def compute_dlog(dlogs, pf, p):
    return sum([dlogs[q] * e for q, e in pf.items()]) % p


def linear_sieve(p, A, B,
                 climit=None, clog=None, eps=None,
                 qlimit=None, fbq=None, fbqlogs=None):
    if climit is None:
        climit = default_climit(p)
    M = A * B // p
    base = A * (B - 1) - M * p
    logbase = base // A + 1
    if clog is None:
        clog = [log(logbase + c) if c > 0 else 0 for c in range(climit)]
    if eps is None:
        eps = 1e-9
    if qlimit is None:
        qlimit = default_qlimit(p)
    if fbq is None:
        fbq = smooth_primes(qlimit)
    if fbqlogs is None:
        fbqlogs = [log(q) for q in fbq]

    sieve = [0] * climit
    num = M * p - A * B
    topp = (M + 1) * p
    topc2 = min(climit, (topp + A - 1) // A - B)
    for q, logq in zip(fbq, fbqlogs):
        if A % q == 0:
            continue
        for yi, qpow in inverse_seq(A, q):
            # assert Hc1 * y % qpow == 1
            c2 = num * yi % qpow
            if A * (B + c2) >= topp or c2 >= climit:
                break
            while c2 < topc2:
                sieve[c2] += logq
                c2 += qpow
    logA = log(A)
    cur = base
    for c2 in range(topc2):
        cur += A
        loglowerbound = logA + clog[c2]
        if loglowerbound < sieve[c2] and abs(sieve[c2] - log(cur)) < eps:
            iss, fn = is_Bsmooth(qlimit, cur)
            if iss:
                yield c2, fn


def worker_linear_sieve(queue, rc2, p, H, clog, eps, qlimit, fbq, fbqlogs):
    for c2 in rc2:
        A, B, climit = H + c2, H, c2 + 1
        res = []
        for c1, fn in linear_sieve(p, A, B, climit,
                                   clog, eps, qlimit, fbq, fbqlogs):
            res.append((c1, fn))
        queue.put((c2, res))
    queue.put(None)


def sign(a):
    if a > 0:
        return 1
    elif a < 0:
        return -1
    else:
        return 0


def ratlift(u, bnd, m):
    a1, a2 = m, u
    v1, v2 = 0, 1
    m2 = bnd
    while True:
        if v2 >= m2:
            return None
        if a2 < m2:
            return sign(v2) * a2, abs(v2)
        q = a1 // a2
        a1 -= q * a2
        v1 -= q * v2
        a1, a2, v1, v2 = a2, a1, v2, v1


def individual_logs(dlogs, Hlogs, y, g, p, op, qlimit, climit, DEBUG=False):
    yy = y
    inf = min(op, 2 ** 31)
    bound = isqrt(p) + 1
    for w in mrange(1, inf, DEBUG=DEBUG):
        yy = yy * g % p
        ab = ratlift(yy, bound, p)
        if ab is not None:
            a, b = ab
            sa = 1 if a > 0 else -1
            a = abs(a)
            iss, af = is_Bsmooth(qlimit, a)
            if iss:
                iss, bf = is_Bsmooth(qlimit, b)
                if iss:
                    ye = (-w + ((p - 1) // 2 if sa == 1 else 0)) % op
                    alog = compute_dlog(dlogs, af, op)
                    blog = compute_dlog(dlogs, bf, op)
                    ye = (ye + alog - blog) % op
                    ye = ye * pow((p - 1) // op, -1, op) % op
                    assert pow(g, ye * (p - 1) // op, p) == y
                    return ye


def individual_logs0(dlogs, Hlogs, y, g, p, op, qlimit, climit, DEBUG=False):
    # this would require ecm factoring to be efficient
    if DEBUG:
        print(f'individual log: y={y} g={g} op={op}')
    assert pow(y, op, p) == 1
    U = L(p, 2 / 3, 3 ** (-1 / 3))
    H = isqrt(p) + 1

    yy = y
    gr = pow(g, (p - 1) // op, p)
    inf = min(op, 2 ** 31)  # tqdm not able to handle bigger numbers
    for w in mrange(1, inf, DEBUG=DEBUG):
        yy = yy * g % p
        fy = factor(yy)
        if max(fy.keys()) < U:
            ye = op - w
            completed = True
            for qq, qe in fy.items():
                if qq < qlimit:
                    ye = (ye + dlogs[qq] * qe) % op
                else:
                    lowu = ceil(H / qq)
                    found = False
                    for u, uf in linear_sieve(p, 1, lowu, climit ** 2):
                        u += lowu
                        for v, nf in linear_sieve(p, u * qq % p, H, climit):
                            if Hlogs[v] is not None:
                                nlog = compute_dlog(dlogs, nf, op)
                                ulog = compute_dlog(dlogs, uf, op)
                                vlog = Hlogs[v]
                                found = True
                                break
                        if found:
                            break
                    if not found:
                        completed = False
                        break
                    ye = (ye + (nlog - ulog - vlog) * qe) % op
            if not completed:
                continue
            # this may not be possible, not sure what could be done
            ye = ye * pow((p - 1) // op, -1, op) % op
            if DEBUG:
                print(f'found log_{gr}({y}) = {ye}')
            assert pow(gr, ye, p) == y
            return ye


def default_climit(p):
    return int(3 * L(p, 0.5, 0.45)) + 10


def default_qlimit(p):
    return int(2 * L(p, 0.5, 0.475)) + 10


class CongruenceFinder:
    def __init__(self, H, p, opq, fbq, fbqlogs, ifbq, qlimit, np):
        self.H, self.p, self.np = H, p, np
        self.fbq, self.fbqlogs, self.qlimit = fbq, fbqlogs, qlimit
        self.opq = opq
        self.c2 = 0
        self.fb = []
        self.infb = []
        self.clog = []
        self.ifbq = ifbq
        self.eps = 1e-9

    def sieve_values(self, n, DEBUG=False):
        for c2 in mrange(self.c2, self.c2 + n, DEBUG=DEBUG):
            for c1, fn in linear_sieve(self.p, self.H + c2, self.H, c2 + 1,
                                       self.clog, self.eps,
                                       self.qlimit, self.fbq, self.fbqlogs):
                yield c2, c1, fn

    def sieve_values_parallel(self, n, DEBUG):
        if (self.c2 + n) <= 1000 or cpu_count() <= 1:
            yield from self.sieve_values(n, DEBUG)
            return
        rng = range(self.c2, self.c2 + n)
        res = parallel_for_balanced(worker_linear_sieve,
                                    (self.p, self.H, self.clog, self.eps,
                                     self.qlimit, self.fbq, self.fbqlogs),
                                    rng, DEBUG=DEBUG)
        for c2 in rng:
            for c1, fn in res[c2]:
                yield c2, c1, fn

    def get(self, n, DEBUG=False):
        if DEBUG:
            print(f'Searching congruences in range {self.c2}-{self.c2 + n}')
        rels, relsex = [], []
        self.infb += [None] * n
        self.clog += [log(c) if c > 0 else 0
                      for c in range(3 * self.c2, 3 * (self.c2 + n))]

        for c2, c1, fn in self.sieve_values_parallel(n, DEBUG):
            if self.infb[c1] is None:
                self.infb[c1] = len(self.fb)
                self.fb += [c1]
            if self.infb[c2] is None:
                self.infb[c2] = len(self.fb)
                self.fb += [c2]
            cr = [0] * self.np
            if c1 == c2:
                crex = [self.infb[c1]]
                for q, e in fn.items():
                    cr[self.ifbq[q]] = e
            else:
                for q, e in fn.items():
                    cr[self.ifbq[q]] = 2 * e
                crex = [self.infb[c1], self.infb[c2]]
            rels += [cr]
            relsex += [crex]
        self.c2 += n
        if DEBUG:
            print(f'resulted in {len(rels)} new rels')
        return rels, relsex


def solve_gaussian(M, q, DEBUG=False):
    row_reduce(M, q, DEBUG)
    m = len(M[0])
    x = [None] * m

    for i in range(m - 1):
        df = 0
        for j in range(m - 1):
            if M[i][j] != 0:
                df += 1
        if df == 1 and M[i][i] == 1:
            x[i] = (-M[i][-1]) % q
    return x


def solve_lanczos(M, q, DEBUG=False):
    b = [(-v[-1]) % q for v in M]
    M, d = to_sp([v[:-1] for v in M])
    x = block_lanczos(M, d, b, q, DEBUG)
    return x


def linear_sieve_dlog(p, gy, y, op=None, qlimit=None, climit=None, DEBUG=False):
    if DEBUG:
        print(f'linear sieve dlog: p={p} gy={gy} y={y}')
    assert isprime(p)
    if op is None:
        op = order(gy, p)
    assert isprime(op)
    assert (p - 1) % op == 0
    assert op >= 10 ** 6
    assert order(gy, p) == op
    opq = op
    while (p - 1) % (opq * op) == 0:
        opq *= op

    if qlimit is None:
        qlimit = default_qlimit(p)
    if climit is None:
        climit = default_climit(p)

    g = primitive_root(p)

    H = isqrt(p) + 1
    assert H > g
    assert H + climit < p

    fbq = smooth_primes(qlimit)
    np = len(fbq)
    ifbq = [None] * qlimit
    for i, q in enumerate(fbq):
        ifbq[q] = i
    fbqlogs = [log(q) for q in fbq]

    rels = []
    relsex = []
    assert g <= qlimit
    if g not in fbq:
        gf = factor(g)
        ig = np
        np += 1
        cr = [0] * np
        for q, e in gf.items():
            assert q < qlimit
            cr[ifbq[q]] = e
        cr[-1] = -1
        rels += [cr]
        relsex += [[]]
    else:
        ig = ifbq[g]

    solved = False
    CF = CongruenceFinder(H, p, opq, fbq, fbqlogs, ifbq, qlimit, np)
    nclimit = None
    first = True
    rounds = 0
    while not solved:
        rounds += 1
        if rounds > 10:
            if DEBUG:
                print('solution not converging, trying higher qlimit')
            return linear_sieve_dlog(p, gy, y, op, qlimit + 50, climit, DEBUG)
        nclimit = climit if first else max(50, climit // 10)
        if DEBUG:
            print(f'\nSolving: p={p} gy={gy} y={y} fbp={len(fbq)}',
                  f'qlimit={qlimit}',
                  f'climit={climit + (0 if nclimit == climit else nclimit)}')
        nrels, nrelsex = CF.get(nclimit, DEBUG)
        rels += nrels
        relsex += nrelsex
        climit += nclimit if not first else 0
        first = False
        if len(nrels) == 0:
            continue
        m = np
        for i in range(climit):
            if CF.infb[i]:
                m += 1
        n = len(rels)
        if n < m + 1:
            continue
        if DEBUG:
            print(f'n={n} relations with m={m} variables')

        rf = structured_gaussian_elimination(rels, relsex, len(CF.fb))
        mod_matrix(rels, opq)
        # rf = rels_filter(rels, relsex, len(CF.fb))
        Mrels, kmap = rels_matrix(rels, relsex, len(CF.fb), rf, DEBUG)
        if len(Mrels) == 0:
            continue
        n, m = len(Mrels), len(Mrels[0]) - np
        if DEBUG:
            print(f'filtered: n={n} relations with m={m + np} variables')
        if len(Mrels) < len(Mrels[0]) + 1:
            continue
        ikmap = [None] * (len(Mrels[0]) - np)
        for i, v in enumerate(kmap):
            if v is not None:
                ikmap[v] = i

        for i in range(n):
            # for j in range(m + np):
            #     Mrels[i][j] %= opq
            # swap g column to the end
            Mrels[i][m + ig], Mrels[i][-1] = Mrels[i][-1], Mrels[i][m + ig]

        # x = solve_gaussian(Mrels, opq, DEBUG)
        x = solve_lanczos(Mrels, opq, DEBUG)

        if x is None:
            continue

        iinfb = [None] * len(CF.fb)
        for i, v in enumerate(CF.infb):
            if v is not None:
                iinfb[v] = i

        gr = pow(g, (p - 1) // opq, p)
        exps = {}
        exps[g] = 1
        Hexps = [None] * climit
        solved = True
        for i in range(m + np - 1):
            if x[i] is not None:
                if i < m:
                    v = iinfb[ikmap[i]]
                    Hexps[v] = x[i]
                elif i - m != ig:
                    exps[fbq[i - m]] = x[i]
                else:
                    exps[fbq[-1]] = x[i]
            elif i >= m:
                # print('more relations required')
                solved = False
                # break

        if solved:
            for i, (rel, relex) in enumerate(zip(rels, relsex)):
                if not rf(i):
                    unknown = [c for c in relex if Hexps[iinfb[c]] is None]
                    if len(unknown) == 1:
                        other = 0
                        for c in relex:
                            if Hexps[iinfb[c]] is not None:
                                other = Hexps[iinfb[c]]
                        v = iinfb[unknown[0]]
                        pf = {fbq[j]: rel[j] for j in range(np) if rel[j] != 0}
                        Hexps[v] = (compute_dlog(exps, pf, opq)
                                    * (opq + 1) // 2 - other) % opq
        if DEBUG:
            if not solved:
                print(f'more relations required {len(exps)}/{np}')

    for q, e in exps.items():
        assert pow(gr, e, p) == pow(q, (p - 1) // opq, p)
    Hlogs = 0
    for v, e in enumerate(Hexps):
        if e is not None:
            Hlogs += 1
            assert pow(gr, e, p) == pow(H + v, (p - 1) // opq, p)

    if DEBUG:
        print(f'{len(exps)}/{np} discrete logs found')
        print(f'H {Hlogs}/{climit} logs found')

    xlogs = []
    for x in [gy, y]:
        xe = individual_logs(exps, Hexps, x, g, p, opq,
                             qlimit, climit, DEBUG=DEBUG)
        if opq != op:
            assert xe % (opq // op) == 0
            xe //= (opq // op)
            assert pow(g, (p - 1) // op * xe, p) == x
        xlogs += [xe]
    ye = xlogs[1] * pow(xlogs[0], -1, op) % op
    assert pow(gy, ye, p) == y
    return ye
