#Use Python3.8
#Gen as subroutine of Setup

import math
import time
import random
import numpy as np
from Crypto.Util import number
from Crypto.Hash import SHAKE256

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

    def Gen(pp,t):
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
    #t = 1000000
    C, t = Gen(pp, t)
    return pp, C, t

#Verifier Encrypt
def Enc(m, pp, e): #e is typically 65537, the most common RSA exponent due to low Hamming weight
    N = pp
    c = pow(m, e, N) #textbook RSA, no OEAP padding https://en.wikipedia.org/wiki/Optimal_asymmetric_encryption_paddingj
    return c

def Enc_OAEP(m, pp, e):                                                         #Completely replaces Enc
    n = 1024                                                                    #Number of bits in the RSA modulus (double the parameter below)
    k0 = 32                                                                     #At least as big as integer type
    k1 = 8
    #m - plaintext message, n - k0 - k1 bits long
    m_prime = (m << k1).to_bytes(int((n-k0)/8), "big")                          #Zero padd message to                n - k0 = 256 bits
    r = random.randint(0, (1<<k0) - 1).to_bytes(int(k0/8), "big")               #Generate a random bit string                  k0 bits
    G = SHAKE256.new(data=r)                                                    #Pass r pass as k0/8 bytes
    Gr = G.read(length=int((n-k0)/8))                                           #Expand r to                         n - k0 = 256 bits
    X = bytes(a ^ b for a, b in zip(m_prime, Gr))                               #XOR the padded message and hash     n - k0 = 256 bits
    H = SHAKE256.new(data=X)                                                    #This reduces 256 bits to k0 bits
    HX = H.read(length=int(k0/8))                                               #Reduce X to                                   k0 bits
    Y = bytes(a ^ b for a, b in zip(r, HX))                                     #XOR r with the new hash                       k0 bits
    X = int.from_bytes(X, "big")
    Y = int.from_bytes(Y, "big")
    msg = (X << k0) | Y                                                         #Concatenate the results
    return Enc(msg, pp, e)                                                      #Encrypt using textbook RSA

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

#Prover Decrypt
def Dec(C, y, c, pp, e):
    N = pp
    p = y[0]
    q = y[1]
    carmichael = np.lcm.reduce([p - 1, q - 1])
    d = pow(e, -1, carmichael) #RSA secret key
    m = pow(c, d, N) #textbook RSA decrypt no OAEP padding          
    return m

def Dec_OAEP(C, y, c, pp, e):
    msg = Dec(C, y, c, pp, e)                                                   #Decrypt using textbook RSA
    n = 1024                                                                    #Modulus size
    k0, k1 = 32, 8                                                              #Integers fixed by protocol
    X = (msg >> k0).to_bytes(int((n-k0)/8), "big")                              #Recover X
    Y = (msg & ((1<<k0)-1)).to_bytes(int(k0/8), "big")                          #Recover Y
    H = SHAKE256.new(data=X)                                                    #Take has of X
    HX = H.read(length=int(k0/8))                                               #to get the value of HX
    r = bytes(a ^ b for a, b in zip(Y, HX))                                     #XOR with Y to recover r
    G = SHAKE256.new(data=r)                                                    #Take the hash of r
    Gr = G.read(length=int((n-k0)/8))                                           #To get the value of Gr
    m_prime = bytes(a ^ b for a, b in zip(X, Gr))                               #XOR with X to recover the padded message
    m_prime = int.from_bytes(m_prime, "big")                                    #Convert back to an integer
    m = m_prime >> k1                                                           #Remove the padding
    return m                                                                    #Return the plaintext message

# run functions
z = 0
t = 100
while z < 4: 
    start_time = time.time()
    pp, C, t = Setup(1,512,t)
    print('\n')
    print('len(N)', len(bin(pp)) - 2)
    #print('N:', pp, '\nC:',  C,  '\nt:', t)
    print('Set:' , round(time.time() - start_time , 4), 'seconds')

    start_time = time.time()
    c = Enc_OAEP(88888, pp, 65537) #change first parameter if you want a different message to be Enc
    #print('c:', c)
    print('Enc:' , round(time.time() - start_time , 4), 'seconds')

    start_time = time.time()
    y = Eval(pp, C, t)
    #y = (1, pp) #lazy eval
    #print('\ny:', y)
    print('Eva:' , round(time.time() - start_time , 4), 'seconds')

    start_time = time.time()
    m = Dec_OAEP(C, y, c, pp, 65537)
    print('m:', m)
    print('Dec:' , round(time.time() - start_time , 4), 'seconds')
    z += 1
