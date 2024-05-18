#! /usr/bin/env python

import pytest
from primefac import isprime
from colorama import Fore

from discretelog.common import factor, random_prime
from discretelog.utils import mrange


def single_test_factor(n):
    print(Fore.GREEN + f"testing n={n}" + Fore.RESET)
    fn = factor(n, DEBUG=True)
    nn = 1
    for p, e in fn.items():
        nn *= p**e
        assert isprime(p)
    assert nn == n


def single_test_factor_random(d, random):
    n = random.randint(10**d, 10 ** (d + 1))
    single_test_factor(n)


def single_test_factor_rsa(d, random):
    n = random_prime(d, random) * random_prime(d, random)
    single_test_factor(n)


def test_factor_small(frandom):
    print("\n" + Fore.RED + "testing factoring algorithms small" + Fore.RESET + "\n")
    for d in mrange(10, 30, 1, True):
        single_test_factor_random(d, frandom)


@pytest.mark.slow
def test_factor_large(frandom):
    print("\n" + Fore.RED + "testing factoring algorithms large" + Fore.RESET + "\n")
    for d in mrange(30, 100, 10, True):
        single_test_factor_random(d, frandom)
    for d in mrange(30, 40, 1, True):
        single_test_factor_rsa(d, frandom)
