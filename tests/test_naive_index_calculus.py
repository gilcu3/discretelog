import pytest

from discretelog.common import random_sophie_germain_prime, primitive_root
from discretelog.utils import mrange
from discretelog.naive_index_calculus import dlog_prime


def single_test(psize, random):
    p = random_sophie_germain_prime(psize, random)
    q = (p - 1) // 2
    b = primitive_root(p) ** 2 % p
    e = random.randint(1, q)
    h = pow(b, e, p)
    assert dlog_prime(b, h, p, True) == e


@pytest.mark.slow
def test_dlog_prime_special():
    tests = [(1540571422742786915303, 25, 690483026481419643586, 641629670911834423534)]
    for p, b, e, h in tests:
        assert dlog_prime(b, h, p, True) == e


def test_dlog_prime_small(frandom):
    ntests = 20
    for psize in [5, 6, 7, 8]:
        for _ in mrange(ntests, DEBUG=True):
            single_test(psize, frandom)


def test_dlog_prime_medium(frandom):
    single_test(15, frandom)
