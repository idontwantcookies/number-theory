from functools import lru_cache
from math import sqrt, floor, isqrt, ceil, log10
from random import getrandbits, randint
from fractions import Fraction


def gcd_extended(a:int, b:int):
	# extended GCD algorith. Time complexity: O(log(min(a,b))
	if a < b:
		r, x, y = gcd_extended(b, a)
		return r, y, x
	r0, r1, r2 = a, b, 0
	x0, x1, x2 = 1, 0, 0
	y0, y1, y2 = 0, 1, 0
	if r1 == 0:
		return r0, x0, y0
	while True:
		r2 = r0 % r1
		q = r0 // r1
		x2 = x0 - q * x1
		y2 = y0 - q * y1
		if r2 == 0:
			return r1, x1, y1
		else:
			r0, x0, y0 = r1, x1, y1
			r1, x1, y1 = r2, x2, y2

def prime(n):
	# deterministic prime in O(sqrt(n))
	n = abs(n)
	if n < 2: return False
	if n < 4: return True
	if n % 2 == 0: return False
	if n < 9: return True
	if n % 3 == 0: return False
	r = isqrt(n)
	for f in range(5, r + 1, 6):
		if n % f == 0: return False
		if n % (f + 2) == 0: return False
	return True

def prime_fermat(n):
	# deterministic Fermat prime test in O(sqrt(n)). faster than prime(), but still slow for big n
	n = abs(n)
	if n == 2: return True
	if n % 2 == 0: return False
	x = isqrt(n)
	if x * x == n: return False
	ceiling = (n + 1) / 2
	while True:
		x += 1
		y = (x * x - n)**0.5
		if x == ceiling: return True
		if y % 1 == 0: return False

def modinv(a, n):
	# returns b such that a * b = 1 (mod n). Time complexity: same as gcd_extended.
	if -2 < n < 2: return 0
	gcd, alfa, _beta = gcd_extended(a, n)
	if gcd != 1:
		return 0
	return alfa % n

def powmod(b, e, n):
	# returns b**e (mod n) using binary exponentiation. Time complexity: O(log(n))
	if abs(n) < 2: raise ValueError(f'n must be an integer with abs(n) > 1.')
	A, P, E = b, 1, e
	while E != 0:
		if E % 2 == 1:
			P = (A * P) % n
		E //= 2
		A = (A * A) % n
	return P

def pre_miller(n:int):
	# encontrando k e q tais que n - 1 = 2^k * q. O(log(n))
	q = n - 1
	k = 0
	while q % 2 == 0:
		q //= 2
		k += 1
	return k, q

def prime_miller(n:int, b:int, k: int, q: int):
	''' Testa se n é um primo usando o teste de miller em base b em O(log(n)). Complexidade: O(k)
	Retorna True caso o número seja TALVEZ primo (inconclusivo), ou 
	False caso o número seja CERTAMENTE composto.
	'''
	if n == 2 or n == -2: return True
	if n % 2 == 0: return False
	if gcd_extended(n, b)[0] == n: return True
	# fazendo o teste de Miller (de fato)
	r = powmod(b, q, n)
	if r == 1 or r == n - 1: return True
	for i in range(1, k):
		r = powmod(r, 2, n)
		if r == n - 1: return True
	return False

def prime_miller_rabin(n:int, rep:int=None):
	''' executa `rep` iterações com bases aleatórioas do teste de Miller 
	para saber se n é primo. Complexidade de tempo: O(rep * log(n))'''
	n = abs(n)
	if n < 2: return False
	if n == 2: return True
	if rep is None:
		rep = ceil(log10(n))
	k, q = pre_miller(n)
	for i in range(rep):
		b = randint(2, n - 1)
		if not prime_miller(n, b, k, q):
			return False
	return True

def random_prime(b:int):
	# retorna um primo aleatório no intervalo [2, 2**b) Complexidade de tempo: O(??).
	# TODO: improve time complexity to not use "bogosort" randomness
	if b < 2: raise ValueError(f'b must be >= 2, not {b}.')
	x = getrandbits(b) | 1 	# force odd number
	while not prime_miller_rabin(x):
		x += 2
	return x

def eratosthenes_sieve(n):
	'''
	Returns a list of all primes between 2 and n. Time complexity: O(n * log(log(n)))
	'''
	l = [True] * max(n, 2)
	out = []
	l[0], l[1] = False, False
	max_i = isqrt(n) + 1
	for i in range(2, max_i):
		if l[i]:
			for j in range(i**2, n, i):
				l[j] = False
	for i, prime in enumerate(l):
		if prime: out.append(i)
	return out

def totient(x, primes=None):
	if prime_miller_rabin(x): return x - 1
	if primes is None: primes = eratosthenes_sieve(x // 2)
	out = x
	for p in primes:
		if x % p == 0:
			out *= (1 - 1 / p)
		if p > x: break
	return floor(out)
