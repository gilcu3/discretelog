#! /usr/bin/env python
import pytest
from discretelog import dlog
from discretelog.utils import mrange
from colorama import Fore


def single_test_random(nb, random):
    n = random.randint(2 ** nb, 2 ** (nb + 1))
    g = random.randint(2, 10)
    e = random.randint(5, n // 2)
    v = pow(g, e, n)
    print(Fore.GREEN + f'testing g={g} e={e} v={v} n={n}' + Fore.RESET)
    d, dt = dlog(g, v, n, DEBUG=True)
    if dt == 0:
        assert d == e, f'd={d}'
    elif dt == 1:
        assert e >= d, f'd={d}'
    else:
        assert e >= d and (e - d) % dt == 0, f'd={d}'


def test_dlog_small(frandom):
    ntests = 1000
    print('\n' + Fore.RED + 'Testing small values' + Fore.RESET + '\n')
    for nb in mrange(5, 30, 1, True):
        for _ in range(ntests):
            single_test_random(nb, frandom)


def test_dlog_medium(frandom):
    print('\n' + Fore.RED + 'Testing medium values' + Fore.RESET + '\n')
    for nb in mrange(30, 60, 1, True):
        single_test_random(nb, frandom)


@pytest.mark.slow
def test_dlog_large(frandom):
    print('\n' + Fore.RED + 'Testing large values' + Fore.RESET + '\n')
    for nb in mrange(80, 200, 5, True):
        single_test_random(nb, frandom)
