#Use Python3.8
#Gen as subroutine of Setup

import math
import time
import random
from Crypto.Util import number

#Verifier Setup - Setup the modulus a Blum Integer
def Setup(x, bits, t):
	def genPrime(bits):
		potential_prime = 1
		while potential_prime % 4 == 1:
			potential_prime = number.getPrime(bits)
		return potential_prime

	x = 0
	while x <1: 
		p = genPrime(bits)
		q = genPrime(bits)
		if p != q and q % 4 != 1:
			N = p * q
			pp = N
		x += 1
	return pp, p , q

#Verifier Gen - even though single use separate Setup and Gen for sequentiality definition
def Gen(pp,t, p, q):
	N = pp
	J_p, J_q = 1, 1
	while not  (J_p == 1 and J_q != 1) and not (J_q == 1 and J_p != 1):
		x = random.randrange(2,N)
		J_p = pow(x,(p-1)//2,p)             #Always == 1 or == p-1, use Euler's Criterion
		J_q = pow(x,(q-1)//2,q)             #Always == 1 or == q-1, use Euler's Criterion
	x_0 = pow(x,2,N)
			
	#now generate x_minus_t
	omega_p = (p + 1) // 4  #Tonelli Shanks, need p = 3 mod 4, Extend Eulers Criterion for proof
	omega_q = (q + 1) // 4  #Tonelli Shanks, need q = 3 mod 4, Extend Eulers Criterion for proof
	
	alpha_p = pow(x_0, pow(omega_p,t,p-1), p) #reduce mod p-1 is Eulers Theorem
	alpha_q = pow(x_0, pow(omega_q,t,q-1), q) #reduce mod p-1 is Eulers Theorem
	
	x_minus_t = ((alpha_p * q * pow(q,-1,p)) + (alpha_q * p * pow(p,-1,q))) % N #Chinese Remainder Theorem to find Mod N. 
	C = (x, x_0, x_minus_t)
	return C, t

#Prover Eval
def Eval(pp, C, t):
	x = C[0]
	x_0 = C[1]
	x_minus_t = C[2]
	N = pp
	
	#evaluation of VDF here
	x_prime = pow(x_minus_t, pow(2,t-1), N) #hard work here
	
	#now use EEA to find the factors of N
	factor1 = math.gcd(x-x_prime, N)
	#factor2 = math.gcd(x+x_prime, N) #O(M(N)logN)
	factor2 = N//factor1 #O(M(N)), so logN better than finding gcd
	y = (factor1, factor2)
	return y

#Verifier Verify
def Verify(pp, C, t, y):
	V = 'reject'
	factor1 = y[0]
	factor2 = y[1]
	x_0 = C[1]
	x_minus_t = C[2]
	N = pp
	
	if factor1 == 1 or factor2 == 1:
		return V
	if (factor1 * factor2 == pp) and (pow(x_minus_t, pow(2,t,(factor1-1)*(factor2-1)) ,N) == x_0):
		V = 'accept'
	return V

# run functions
z = 0
t = 1000000
#t=10
while z < 2: 
	start_time = time.time()
	pp, p , q = Setup(1,1024,t)
	print('\npp:', pp)
	print('Set:' , round(time.time() - start_time , 4), 'seconds')
	
	start_time = time.time()
	C, t = Gen(pp,t, p, q)
	print('C, t:', C, t)
	print('Gen:' , round(time.time() - start_time , 4), 'seconds')

	start_time = time.time()
	y = Eval(pp, C, t)
	#y = (1, pp) #lazy eval
	print('y:', y)
	print('Eva:' , round(time.time() - start_time , 4), 'seconds')

	start_time = time.time()
	V = Verify(pp, C, t, y)
	print('V:', V)
	print('Ver:' , round(time.time() - start_time , 4), 'seconds')
	z += 1
