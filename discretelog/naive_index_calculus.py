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
from .utils import mrange, NoOpTqdm


def producttree(X: list[int]) -> list[list[int]]:
    result = [X]
    while len(X) > 1:
        X = [prod(X[i * 2 : (i + 1) * 2]) for i in range((len(X) + 1) // 2)]
        result.append(X)
    return result


def remaindersusingproducttree(n: int, T: list[list[int]]) -> list[int]:
    result = [n]
    for t in reversed(T):
        result = [result[i // 2] % t[i] for i in range(len(t))]
    return result


def remainders(n: int, X: list[int]) -> list[int]:
    return remaindersusingproducttree(n, producttree(X))


def berstein_prec(b: int) -> int:
    ps = smooth_primes(b)
    z = prod(ps)
    return z


congruence = tuple[dict[int, int], int]


# this is probably not worth it in python
def berstein_batch_smooth(
    b: int, ns: list[tuple[int, int]], z: int
) -> list[congruence]:
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
    def __init__(self, g: int, B: int, p: int):
        self.g, self.B, self.p = g, B, p
        self.congs: list[congruence] = []
        self.k = 0

    def get(self, n: int, DEBUG: bool = False) -> tuple[list[int], list[congruence]]:
        if DEBUG:
            print(f"Searching {n} congruences")
        batch_size = self.B * 4
        z = berstein_prec(self.B)

        pbar = tqdm(total=n) if DEBUG else NoOpTqdm()

        found = 0
        while found < n:
            bk = []
            for _ in range(batch_size):
                self.k += 1
                bk += [(pow(self.g, self.k, self.p), self.k)]
            pss = berstein_batch_smooth(self.B, bk, z)
            found += len(pss)
            self.congs += pss

            pbar.update(len(pss))
        pbar.close()
        bases = list(set(base for c in self.congs for base in c[0]))
        if DEBUG:
            print(f"resulted in bases: {len(bases)}")
        return bases, self.congs


def to_matrices(bases: list[int], congruences: list[congruence]) -> list[list[int]]:
    M = [
        [c[0][base] if base in c[0] else 0 for base in bases] + [c[1]]
        for c in congruences
    ]
    return M


def modMatrix(xx: list[list[int]], m: int) -> list[list[int]]:
    return [[elem % m for elem in sublist] for sublist in xx]


def evaluate(eq: dict[int, int], dlogs: dict[int, int], q: int) -> int:
    return sum([dlogs[term] * exp for term, exp in eq.items()]) % q


def check_congruences(
    congruences: list[congruence], dlogs: dict[int, int], q: int
) -> bool:
    passed = True
    for c in congruences:
        if evaluate(c[0], dlogs, q) != c[1] % q:
            passed = False
            break
    return passed


def check_dlogs(g: int, p: int, q: int, exponents: list[int], bases: list[int]) -> bool:
    passed = True
    o = phi(p)
    for exponent, base in zip(exponents, bases):
        if pow(g, exponent * o // q, p) != pow(base, o // q, p):
            passed = False
            break
    return passed


def msolve_prime(M: list[list[int]], q: int, DEBUG: bool = False) -> list[int] | None:
    if DEBUG:
        print(f"solving linear system {len(M)}x{len(M[0])}")
    n = len(M[0]) - 1
    m = modMatrix(M, q)
    if not row_reduce(m, q):
        return None
    return [m[i][n] for i in range(n)]


def dlog_prime(b: int, h: int, p: int, DEBUG: bool = False) -> int:
    if p <= 10**5:
        raise AssertionError(f"{p} is too small")
    q = multiplicative_order(b, p)
    assert isprime(q)
    o = phi(p)
    g = primitive_root(p)
    B = int(4 * exp(sqrt(log(p) * log(log(p))) / 2)) + 10
    if DEBUG:
        print("p: {}, b: {}, h: {}, B: {}".format(p, b, h, B))
    congs: list[congruence] = []
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

    def find_log_g(x: int) -> int:
        if DEBUG:
            print(f"searching for k such that {x}*g^-k is B-smooth.")
        for k in mrange(1, 10**7, 1, DEBUG):
            c = is_Bsmooth(B, x * pow(g, -k, p) % p)
            if c[0] and set(c[1].keys()).issubset(set(dlogs.keys())):
                if DEBUG:
                    print("found k = {}".format(k))

                xe = (evaluate(c[1], dlogs, q) + k) % q
                if DEBUG:
                    print(f"Found partial dlog: {g}^{xe}={x} mod {p} mod {q}")
                assert pow(g, xe * o // q, p) == pow(x, o // q, p)
                return xe
        raise Exception(f"find_log_g could not find solution {x=} {g=} {p=}")

    be = find_log_g(b)
    he = find_log_g(h)

    sol = he * pow(be, -1, q) % q
    if pow(b, sol, p) == h:
        if DEBUG:
            print("{}^{} = {} (mod {}) holds!".format(b, sol, h, p))
            print("DLP solution: {}".format(sol))
        return sol

    raise Exception(f"dlog_prime could not find solution {b=} {h=} {p=}")
