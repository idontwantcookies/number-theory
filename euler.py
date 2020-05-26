'''
I developed most of these while studying number theory and/or solving problems
on the subject. Many of these functions already have a python built-in solution,
but I wrote raw code anyway just to see if I could. Other functions, however, 
are legitimately useful and are not found in python's base lib, although I have 
found third-party implementations whose focus was on number theory.
'''

from functools import lru_cache
from math import sqrt, floor, isqrt, ceil, log10
from random import getrandbits, randint

def digits(n):
	n = abs(n)
	out = []
	while n > 0:
		out.insert(0, n % 10)
		n //= 10
	return out

def digit_sum(x):
	out = 0
	while x > 0:
		out += x % 10
		x //= 10
	return out

def reverse(n):
	out = 0
	while n > 0:
		out = 10 * out + n % 10
		n //= 10
	return out

def sum_multiples(n, target):
	p = target // n
	return n * (p * (p + 1) // 2)

def sum_naturals(n):
	# returns the sum of n first natural numbers
	return n * (n + 1) // 2

def sumsq_naturals(n):
	# returns the sum of the squares of the first n natural numbers
	return (2 * n + 1) * (n + 1) * n // 6

@lru_cache(maxsize=None)
def fibonacci(n):
	if n == 0: return 0
	if n == 1: return 1
	return fibonacci(n - 1) + fibonacci(n - 2)

def divide_until_not_multiple(n, factor):
	while n % factor == 0:
		n //= factor
	return n

def largest_prime_factor(n):
	n = divide_until_not_multiple(n, 2)
	max_factor = sqrt(n)
	factor = 3
	while n > 1 and factor <= max_factor:
		if n % factor == 0:
			n = divide_until_not_multiple(n, factor)
			last_factor = factor
			max_factor = sqrt(n)
		factor += 2
	return last_factor if n == 1 else n

def is_palindrome(n):
	return n == reverse(n)

def prime(n):
	# deterministic prime (very slow for big numbers)
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
	# deterministic Fermat prime test (faster than prime(), but still slow for big n)
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
	# returns b such that a * b = 1 (mod n)
	if -2 < n < 2: raise ValueError(f'n must be an integer greater than 1.')
	gcd, alfa, beta = gcd_extended(a, n)
	if gcd != 1:
		raise ValueError(f'{a} does not have a modular inverse in {n} because gcd({a},{n}) = {gcd} != 1.')
	return alfa % n

def powmod(b, e, n):
	# returns b**e (mod n) using binary exponentiation.
	if abs(n) < 2: raise ValueError(f'n must be an integer with abs(n) > 1.')
	A, P, E = b, 1, e
	while E != 0:
		if E % 2 == 1:
			P = (A * P) % n
		E //= 2
		A = (A * A) % n
	return P

def pre_miller(n:int):
	# encontrando k e q tais que n - 1 = 2**k * q
	n1 = n - 1
	q = n1
	k = 0
	while True:
		q //= 2
		k += 1
		if q % 2 == 1: break
	return n1, k, q

def prime_miller(n:int, b:int, n1: int, k: int, q: int):
	''' testa se n é um primo usando o teste de miller em base b. retorna
	True caso o número seja TALVEZ primo (inconclusivo), ou 
	False caso o número seja CERTAMENTE composto.'''
	if n == 2 or n == -2: return True
	if n % 2 == 0: return False
	if gcd(n, b) == n: return True
	# fazendo o teste de Miller (de fato)
	r = powmod(b, q, n)
	if r == 1 or r == n1: return True
	for i in range(1, k):
		r = powmod(r, 2, n)
		if r == n1: return True
	return False

@lru_cache(maxsize=None)
def prime_miller_rabin(n:int, iter:int=None):
	''' executa (iter) iterações com bases aleatórioas do teste de Miller 
	para saber se n é primo.'''
	n = abs(n)
	if n < 2: return False
	if n == 2: return True
	if iter is None:
		iter = ceil(log10(n))
	n1, k, q = pre_miller(n)
	for i in range(iter):
		b = randint(2, n1)
		if not prime_miller(n, b, n1, k, q):
			return False
	return True

def random_prime(b:int):
	# retorna um primo aleatório no intervalo [2, 2**b).
	if b < 2: raise ValueError(f'b must be >= 2, not {b}.')
	while True:
		x = getrandbits(b)
		x |= 1	# making sure its an odd number
		if prime_miller_rabin(x): return x

def eratosthenes_siege(n):
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

def polygonal_number(s, n):
	return ((s - 2) * n - (s - 4)) * n // 2

def cfrac_sqrt(x):
	period = []
	m, d = 0, 1
	a0 = isqrt(x)
	if a0 * a0 == x: return a0, []
	a = a0
	while a != 2 * a0:
		m = d * a - m
		d = (x - m * m) / d
		a = int((a0 + m) / d)
		period.append(a)
	return a0, period

def convergent_fractions(a0, period):
	h_1, h_2 = 1, 0
	k_1, k_2 = 0, 1
	a = a0
	i = -1
	N = len(period)
	while True:
		h = a * h_1 + h_2
		k = a * k_1 + k_2
		yield h, k
		h_2, h_1 = h_1, h
		k_2, k_1 = k_1, k
		i += 1
		a = period[i % N]

def e_convergent(depth):
	frac = Fraction(0)
	for i in range(depth, 0, -1):
		mod = i % 3
		k = (i + 1) // 3
		x = 2 * k if mod == 2 else 1
		frac += x
		frac = 1 / frac
	return frac + 2

def partitions(target, array=None):
	target += 1
	if array is None: array = tuple(range(1, target))
	ways = [0] * (target)
	ways[0] = 1
	for x in array:
		for j in range(x, target):
			ways[j] += ways[j - x]
	return ways[-1]

def gcd(a, b):
	if b > a: a, b = b, a
	while b != 0:
		a, b = b, a % b
	return a

def lcm(a, b):
	return a * b // gcd(a, b)

def gcd_extended(a:int, b:int):
	swap = a < b
	r = [a, b, 0] if not swap else [b, a, 0]	# remainders
	alfa = [1, 0, 0]
	beta = [0, 1, 0]
	if r[1] == 0:
		return (r[0], alfa[0], beta[0]) if not swap else (r[0], beta[0], alfa[0])
	while True:
		r[2] = r[0] % r[1]
		quociente = r[0] // r[1]
		alfa[2] = alfa[0] - quociente * alfa[1]
		beta[2] = beta[0] - quociente * beta[1]
		if r[2] == 0:
			return (r[1], alfa[1], beta[1]) if not swap else (r[1], beta[1], alfa[1])
		else:
			for i in range(2):
				(r[i], alfa[i], beta[i]) = (r[i+1], alfa[i+1], beta[i+1])

def totient(x, primes=None):
	if prime_miller_rabin(x): return x - 1
	if primes is None: primes = eratosthenes_siege(x)
	out = x
	for p in primes:
		if x % p == 0:
			out *= (1 - 1 / p)
		if p > x: break
	return int(out)

def pytagorean_triplet(a, b):
	if a < b: a, b = b, a
	x = a**2 - b**2
	y = 2 * a * b
	d = gcd(x, y)
	x, y = x // d, y // d
	z = (a**2 + b**2) // d
	if x > y: x, y = y, x
	return x, y, z
