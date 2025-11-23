import pytest

from discretelog.block_lanczos import spmul2, spmul, spmax
from discretelog.common import random_prime


def random_sparse_matrix(n, s, bound, frandom):
    mat = [[] for _ in range(n)]
    for i in range(n):
        for j in sorted(frandom.choices(range(n), k=s)):
            mat[i].append((j, frandom.randint(-bound, bound)))
    return mat


def random_vector(n, bound, p, frandom):
    return [frandom.randint(-bound, bound) % p for _ in range(n)]


@pytest.mark.perf
def test_spmul2_15(frandom, benchmark):
    n, s, bound = 1000, 100, 30
    p = random_prime(15, frandom)
    mat = random_sparse_matrix(n, s, bound, frandom)
    a = random_vector(n, bound, p, frandom)
    mx = spmax(mat)
    benchmark(spmul2, mx, mat, a, p)


@pytest.mark.perf
def test_spmul_15(frandom, benchmark):
    n, s, bound = 1000, 100, 30
    p = random_prime(15, frandom)
    mat = random_sparse_matrix(n, s, bound, frandom)
    a = random_vector(n, bound, p, frandom)
    benchmark(spmul, mat, a, p)


@pytest.mark.perf
def test_spmul2_16(frandom, benchmark):
    n, s, bound = 1000, 100, 30
    p = random_prime(16, frandom)
    mat = random_sparse_matrix(n, s, bound, frandom)
    a = random_vector(n, bound, p, frandom)
    mx = spmax(mat)
    benchmark(spmul2, mx, mat, a, p)


@pytest.mark.perf
def test_spmul_16(frandom, benchmark):
    n, s, bound = 1000, 100, 30
    p = random_prime(16, frandom)
    mat = random_sparse_matrix(n, s, bound, frandom)
    a = random_vector(n, bound, p, frandom)
    benchmark(spmul, mat, a, p)


@pytest.mark.perf
def test_spmul2_18(frandom, benchmark):
    n, s, bound = 1000, 100, 30
    p = random_prime(18, frandom)
    mat = random_sparse_matrix(n, s, bound, frandom)
    a = random_vector(n, bound, p, frandom)
    mx = spmax(mat)
    benchmark(spmul2, mx, mat, a, p)


@pytest.mark.perf
def test_spmul_18(frandom, benchmark):
    n, s, bound = 1000, 100, 30
    p = random_prime(18, frandom)
    mat = random_sparse_matrix(n, s, bound, frandom)
    a = random_vector(n, bound, p, frandom)
    benchmark(spmul, mat, a, p)
