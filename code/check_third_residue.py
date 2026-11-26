"""Exact verification of the third connected residue.

Run beside verify.py (imported from this directory) with Python 3.10+ and SymPy.
The character-residue algebra is verified symbolically for general d, f, c.
The independent cut-and-join counts verify all profiles of degrees five to eight.
"""
from math import factorial
import sympy as sp
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify import cut_join, fit, integer_partitions, polynomial_A

def check(condition, message=''):
    """Explicit check that also runs under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))
d, f, c, t = sp.symbols('d f c t')
D, u = (d - 1, f * (f - 1))
FD = 1 - sp.Rational(1, 2) * (d - 1) * t
FD1 = 1 - (d - 2) * t
Fminus1 = 1 - t / (d - 2)
Fminus2 = 1 - 2 * t / (d - 1)
Fone = 1 + t / (d - 4)
shape_terms = FD * FD1 + D * (f - 1) * FD1 * Fminus1 + d * (d - 3) / 2 * (f * (f - 1) / 2 + c - f) * Fminus1 + (d - 1) * (d - 2) / 2 * ((f - 1) * (f - 2) / 2 - c) * Fminus1 * Fminus2
singleton_term = -d * f * (FD1 + (d - 2) * (f - 2) * Fminus1)
pair_term = -d * D / 4 * ((u + 2 * c) * Fone + (u - 2 * c) * Fminus1)
two_singleton_term = d * D * u
a = (d + f - 2) * (d + f - 3) / 2 - c
Q = (d - 1) * (2 * c * (d - 2) + (d - 3) * (d - 4) + 2 * f * (d + f - 5)) / ((d - 4) * (d - 2))
check(sp.factor(shape_terms + singleton_term + pair_term + two_singleton_term - a * (1 + t) ** 2 + Q * t) == 0, 'sp.factor(shape_terms + singleton_term + pair_term + two_singleton_term - a * (1 + t) ** 2 + Q * t) == 0')
print('General rational identity in d, f, c, t verified exactly')
g = sp.Symbol('g')
for degree in range(5, 9):
    bases = [(k, degree // (k + 1)) for k in range(1, degree)]
    number = sum((order for k, order in bases))
    profiles = integer_partitions(degree)
    for profile in profiles:
        values = [cut_join(genus, profile) for genus in range(number + 2)]
        r0 = degree + len(profile) - 2
        pole = degree - 3
        spectrum = [fit([v[h] for v in values], bases, g) for h in range(1, degree)]
        residue = sp.expand(sp.Rational(factorial(degree), 2 * pole ** r0) * sum((s[pole * pole] * t ** h for h, s in enumerate(spectrum, 1))))
        substitution = {d: degree, f: profile.count(1), c: profile.count(2)}
        predicted = polynomial_A(pole, t) * (a * (1 + t) ** 2 - Q * t).subs(substitution)
        check(sp.expand(residue - predicted) == 0, (degree, profile))
    print(f'd={degree}: all {len(profiles)} profiles agree with the third residue')
