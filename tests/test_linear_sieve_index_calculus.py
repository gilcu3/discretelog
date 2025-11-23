import pytest
from primefac import isprime
from colorama import Fore

from discretelog.common import (
    random_prime,
    random_sophie_germain_prime,
    primitive_root,
    order,
)
from discretelog.utils import mrange
from discretelog.linear_sieve_index_calculus import linear_sieve_dlog


def single_test(psize, frandom):
    p = random_sophie_germain_prime(psize, frandom)
    q = (p - 1) // 2
    g = primitive_root(p)
    gr = g**2 % p
    e = frandom.randint(2, q)
    y = pow(gr, e, p)
    assert pow(gr, q, p) == 1
    assert pow(y, q, p) == 1
    ye = linear_sieve_dlog(p, gr, y, q, DEBUG=True)
    assert e == ye


def higher_order_test(psize, frandom):
    p = random_prime(psize, frandom)
    while True:
        r = frandom.randint(1, 10**psize)
        if isprime(r * p**2 + 1):
            break
    p, q = r * p**2 + 1, p
    g = primitive_root(p)
    gr = pow(g, q * r, p)
    e = frandom.randint(1, q)
    y = pow(gr, e, p)
    ye = linear_sieve_dlog(p, gr, y, q, DEBUG=True)
    assert e == ye


@pytest.mark.slow
def test_special(frandom):
    for p, g, e, op in [
        (1000000000000000003, None, None, 52445056723),
        (100000000000000000039, None, None, 507526619771207),
        (1000000000000000000000000000057, None, None, 454197539),
        # very slow test for relation finding phase
        (
            33380411190168454492618748210515919,
            28236217495251472904421868195362610,
            1040261549503090216732095527114575838939980980,
            None,
        ),
        # very slow for individual logs phase
        (
            965646041566206382000881149016637,
            338297506802314145874810655853096,
            157024015806831411,
            608402596429614983,
        ),
    ]:
        if g is None or op is None:
            if g is None:
                g = primitive_root(p)

            if op is None:
                op = order(g, p)

            gr = pow(g, (p - 1) // op, p)
        else:
            gr = g
        assert isprime(op)

        if e is None:
            e = frandom.randint(1, op)
        else:
            e %= op

        y = pow(gr, e, p)
        print(Fore.GREEN + f"testing {gr=} {e=} {y=} {p=}" + Fore.RESET)
        ye = linear_sieve_dlog(p, gr, y, op, DEBUG=True)
        assert ye == e


def test_linear_sieve_dlog_small(frandom):
    ntests = 100
    print("sophie germaine primes")
    for psize in range(6, 9):
        print(f"psize={psize}")
        for _ in mrange(ntests, DEBUG=True):
            single_test(psize, frandom)


@pytest.mark.slow
def test_linear_sieve_dlog_medium(frandom):
    ntests = 100
    print("sophie germaine primes")
    for psize in range(9, 14):
        print(f"psize={psize}")
        for _ in mrange(ntests, DEBUG=True):
            single_test(psize, frandom)
    ntests = 100
    print("higher order")
    for _ in mrange(ntests, DEBUG=True):
        higher_order_test(6, frandom)


@pytest.mark.slow
def test_linear_sieve_dlog_large(frandom):
    ntests = 1
    print("sophie germaine primes")
    for psize in range(27, 30):
        print(f"psize={psize}")
        for _ in mrange(ntests, DEBUG=True):
            single_test(psize, frandom)
    ntests = 1
    print("higher order")
    for _ in mrange(ntests, DEBUG=True):
        higher_order_test(10, frandom)
