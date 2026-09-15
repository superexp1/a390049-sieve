#!/usr/bin/env python3
"""Independent reference for A390049, written without reference to the C++.

Factors every n < N with a smallest-prime-factor sieve and evaluates
sigma(n) == psi(n) + phi(n) + omega(n) directly from the factorisation.
Deliberately naive: its only job is to disagree with the sieve if the
sieve is wrong.

    tests/reference.py [N]        # default 3000000
"""
import sys

def terms_below(N):
    spf = list(range(N))
    i = 2
    while i * i < N:
        if spf[i] == i:
            for j in range(i * i, N, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    out = []
    for n in range(2, N):
        m, sigma, psi, phi, omega = n, 1, 1, 1, 0
        while m > 1:
            p, e = spf[m], 0
            while m % p == 0:
                m //= p
                e += 1
            sigma *= (p ** (e + 1) - 1) // (p - 1)
            psi   *= p ** (e - 1) * (p + 1)
            phi   *= p ** (e - 1) * (p - 1)
            omega += 1
        if sigma == psi + phi + omega:
            out.append(n)
    return out

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000000
    for t in terms_below(N):
        print(t)
