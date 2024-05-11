#! /usr/bin/env python

import pytest
from primefac import isprime
from colorama import Fore

from discretelog.common import factor
from discretelog.utils import mrange


def single_test_factor(d, random):
    n = random.randint(2 ** d, 2 ** (d + 1))
    print(Fore.GREEN + f'testing n={n}' + Fore.RESET)
    fn = factor(n, DEBUG=True)
    nn = 1
    for p, e in fn.items():
        nn *= p ** e
        assert isprime(p)
    assert nn == n


def test_factor_small(frandom):
    print('\n' + Fore.RED + 'testing factoring algorithms small'
          + Fore.RESET + '\n')
    for d in mrange(20, 50, 1, True):
        single_test_factor(d, frandom)


@pytest.mark.slow
def test_factor_large(frandom):
    print('\n' + Fore.RED + 'testing factoring algorithms large'
          + Fore.RESET + '\n')
    for d in mrange(50, 350, 10, True):
        single_test_factor(d, frandom)
