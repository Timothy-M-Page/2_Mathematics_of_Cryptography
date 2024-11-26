import math

# To find the greatest common divisor of two integers and inverses
# in Z/mZ, the extended Euclidean algorithm may be implemented:


def extended_euclid(a, b):       # To find the gcd of the integers a and b.
    if a == 0:
        return [b, [0, 1]]
    if b == 0:
        return [a, [1, 0]]       # Handle the case that one or both of the inputs is zero.
    A = max(a, b)
    B = min(a, b)
    r = A % B
    q = A//B
    R = [A, B, r]   # Remainders
    Q = [q]         # Quotients
    S = [1, 0]      # The sets S and T allow the integer coefficients to be calculated iteratively
    T = [0, 1]      # that allow the gcd to be expressed as ax + by = gcd.
    while R[-1] != 0:
        q = R[-2]//R[-1]
        Q.append(q)
        r = R[-2] % R[-1]
        R.append(r)
        s = S[-2] - Q[-2]*S[-1]
        S.append(s)
        t = T[-2] - Q[-2]*T[-1]
        T.append(t)
    if a == A:
        return [R[-2], S[-1], T[-1]]    # Returns the gcd and the integer coefficients allowing the gcd
    else:                               # to be expressed as a linear combination of a and b.
        return [R[-2], T[-1], S[-1]]


def gcd(a, b):
    return extended_euclid(a,b)[0]

# Python implements the extended Euclidean algorithm with the pow function:


def inv(a, mod):
    if gcd(a, mod) > 1:
        return str(a) + " has no inverse mod " + str(mod)
    return pow(a, -1, mod)


# A function essential for RSA is Euler's phi function:

def phi(n):
    result = n         # Start with the value n
    p = 2
    while p**2 <= n:
        if n % p == 0:
            while n % p == 0:
                n //= p
            result -= result // p     # For each prime that divides n, multiply result by (1-1/p)
        p += 1                        # subtracting result//p is equivalent to this multiplication.
    if n > 1:       # A condition for the case that n itself is prime.
        result -= result // n
    return result

# The phi function will later be used in RSA through the application of
# Euler's theorem stating that a**phi(n) = 1 mod(n) for all a such that gcd(a,n) = 1.


# Pollard's Rho Algorithm :

# Cracking RSA involves calculating phi(n) when n is the product of two large primes.
# Factoring numbers into primes is an NP problem, however there are more efficient algorithms
# than simply checking each prime iteratively. Below is an implementation of Pollard's algorithm.

def pollard(n):
    m = n   # Store the original input
    A = []  # An array to collect factors
    x = 2
    y = 2
    i = 0
    def f(x):
        return x**2     # Choose a polynomial
    while n > 1:
        i += 1
        x = f(x) % n
        y = f(f(y)) % n
        d = gcd(abs(x-y), n)    # Iterate the values of x,y and d until a factor is shared with n
        if d == 1:
            continue
        if d in range(2, n):
            n = n//d            # Remove each factor when it is found and continue the loop
            A.append(d)
            if n == 1:          # If all factors are found, break the loop
                break
        if d == n:
            if math.prod(A) == m or math.prod(A)*d == m:  # A condition for if the remaining factor n//d is itself prime
                A.append(d)
                break
            else:
                return "Polynomial failure due to gcd(|x-y|,n) = n at step " + str(i) + "."
    return A



# Chinese remainder theorem:

# For a set of moduli M and a set of remainders R, if P is the product of
# the elements of M. Then there exists one x mod P such that x = R[i] % M[i] for all i.


def chinese_remainder_theorem(moduli, remainders):      # A function of two arrays
    P = math.prod(moduli)
    N = []
    for m in moduli:
        N.append(P//m)                      # Construct the set of products of all moduli but one for each modulus
    I = []
    for i in range(len(moduli)):
        I.append(inv(N[i],moduli[i]))       # Calculate the modular inverses of these values
    x = 0
    for i in range(len(moduli)):
        x = x + (remainders[i]*N[i]*I[i])   # Taking this expression mod moduli[i] necessarily returns remainder[i] by construction
    x = x % P
    return x


