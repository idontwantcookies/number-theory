import pytest

import tp

extgcd = []
with open("test/extended_gcd.txt") as file:
    for line in file:
        extgcd.append([int(x) for x in line.split()])

fastexp = []
with open("test/exp_binaria.txt") as file:
    for line in file:
        fastexp.append([int(x) for x in line.split()])

invmod = []
with open("test/inverso_modular.txt") as file:
    for line in file:
        invmod.append([int(x) for x in line.split()])

primes = []
with open("test/primes.txt") as file:
    for line in file:
        primes.append([int(x) for x in line.split()])

@pytest.mark.parametrize("n,r", [
    [1, 1],
    [2, 1],
    [3, 1],
    [4, 2],
    [5, 2],
    [6, 2],
    [7, 2],
    [8, 2],
    [9, 3],
    [101, 10]
])
def test_isqrt(n, r):
    assert tp.isqrt(n) == r

@pytest.mark.parametrize("a,b,x,y,d", extgcd)
def test_gcd_extended(a, b, x, y, d):
    assert tp.gcd_extended(a, b) == (d, x, y)


@pytest.mark.parametrize("b,e,n,x", fastexp)
def test_exp_binaria(b, e, n, x):
    assert tp.powmod(b, e, n) == x


@pytest.mark.parametrize("a,n,inv", invmod)
def test_inverso_modular(a, n, inv):
    assert tp.modinv(a, n) == inv


@pytest.mark.parametrize("p,result", primes)
def test_primes_miller_rabin(p, result):
    result = bool(result)
    assert tp.prime_miller_rabin(p) == result
