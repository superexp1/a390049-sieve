#!/usr/bin/env python3
"""Independent reference for oddpart, written from the mathematics only.

For odd m, with P,Q,R = sigma,psi,phi of m and o = omega(m):

    odd k = m        is a term iff  P == Q + R + o
    even k = 2^a*m   is a term iff  X = (P + o + 1) / (4P - 3Q - R)
                                    is an exact power of two, X = 2^(a-1)

Emits exactly the lines ./oddpart prints, so the two can be diffed.

    tests/oddpart-reference.py [M]        # default 200000
"""
import sys

def main(M):
    if M < 3:
        return []
    spf = list(range(M))
    i = 2
    while i * i < M:
        if spf[i] == i:
            for j in range(i * i, M, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    out = []
    for m in range(3, M, 2):
        n, P, Q, R, o = m, 1, 1, 1, 0
        while n > 1:
            p, e = spf[n], 0
            while n % p == 0:
                n //= p
                e += 1
            P *= (p ** (e + 1) - 1) // (p - 1)
            Q *= p ** (e - 1) * (p + 1)
            R *= p ** (e - 1) * (p - 1)
            o += 1
        if P == Q + R + o:
            out.append("ODD-TERM k=%d omega=%d" % (m, o))
        D = 4 * P - 3 * Q - R
        if D > 0:
            num = P + o + 1
            if num % D == 0:
                X = num // D
                if X >= 1 and (X & (X - 1)) == 0:
                    out.append("TERM oddpart=%d X=%d omega=%d" % (m, X, o + 1))
    return out

if __name__ == "__main__":
    for line in main(int(sys.argv[1]) if len(sys.argv) > 1 else 200000):
        print(line)
