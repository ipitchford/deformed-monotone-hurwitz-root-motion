"""Exact monotone-transposition word counts by permutation dynamic programming.

Requires Python 3.10+ and NumPy. SymPy isolates real roots to rational intervals.
No characters, contents or generating-function residues are used for the counts.
"""
from itertools import permutations, product
from math import factorial
import numpy as np
import sympy as sp

def check(condition, message=''):
    """Explicit check that also runs under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))

def set_partitions(items):
    if not items:
        yield ()
        return
    first, *rest = items
    for tail in set_partitions(rest):
        yield (((first,),) + tail)
        for k, block in enumerate(tail):
            yield (tail[:k] + ((first,) + block,) + tail[k + 1:])

def restricted_counts(blocks, target, max_length, digit_bits):
    states = []
    for images in product(*(permutations(block) for block in blocks)):
        permutation = list(range(len(target)))
        for block, image in zip(blocks, images):
            for vertex, value in zip(block, image):
                permutation[vertex] = value
        states.append(tuple(permutation))
    index = {p: i for i, p in enumerate(states)}
    nstates = len(states)
    counts = np.zeros((max_length + 1, nstates), dtype=object)
    counts[0, index[tuple(range(len(target)))]] = 1
    for b in range(1, len(target)):
        block = next((block for block in blocks if b in block))
        small = [a for a in block if a < b]
        if not small:
            continue
        exchanges = []
        for a in small:
            predecessors = []
            for perm in states:
                pred = list(perm)
                pred[a], pred[b] = (pred[b], pred[a])
                predecessors.append(index[tuple(pred)])
            exchanges.append(np.asarray(predecessors))
        updated = counts.copy()
        ending = np.zeros(nstates, dtype=object)
        for length in range(1, max_length + 1):
            candidates = ending + counts[length - 1]
            ending = sum((candidates[pred] for pred in exchanges), np.zeros(nstates, dtype=object))
            updated[length] += ending << digit_bits
        counts = updated
    return counts[:, index[target]]

def connected_counts(profile, max_length):
    d = sum(profile)
    cycles, start = ([], 0)
    target = list(range(d))
    for length in profile:
        cycle = tuple(range(start, start + length))
        for i, vertex in enumerate(cycle):
            target[vertex] = cycle[(i + 1) % length]
        cycles.append(cycle)
        start += length
    bits = 1 + max_length * (d * (d - 1) // 2).bit_length()
    packed = np.zeros(max_length + 1, dtype=object)
    for partition in set_partitions(list(range(len(cycles)))):
        blocks = [tuple(sorted((v for i in part for v in cycles[i]))) for part in partition]
        k = len(blocks)
        packed += (-1) ** (k - 1) * factorial(k - 1) * restricted_counts(blocks, tuple(target), max_length, bits)
    mask = (1 << bits) - 1
    return [[p >> bits * i & mask for i in range(d)] for p in packed]
t = sp.Symbol('t')

def leading_poly(k):
    p = t
    for j in range(1, k):
        p *= 1 + sp.Rational(j, k - j) * t
    return sp.Poly(p, t)

def check_profile(profile, max_length):
    d = sum(profile)
    D, E = (d - 1, d - 2)
    m = profile.count(1)
    counts = connected_counts(profile, max_length)
    residuals = []
    parity = (d - len(profile)) % 2
    first = 8 + parity
    for r in range(first, max_length + 1, 2):
        direct = sp.Poly(sum((c * t ** j for j, c in enumerate(counts[r]))), t)
        predicted = D ** r * leading_poly(D) - (E + m) * E ** r * sp.Poly(1 + t, t) * leading_poly(E)
        residuals.append(factorial(d) * direct / 2 - predicted.as_expr())
    for j in range(1, E):
        for power in range(D + 1 - j):
            residuals = [sp.expand(b - j * j * a) for a, b in zip(residuals, residuals[1:])]
    check(len(residuals) >= 2, 'len(residuals) >= 2')
    check(all((sp.expand(q) == 0 for q in residuals)), 'all((sp.expand(q) == 0 for q in residuals))')
    print(f'{profile}: {len(residuals)} exact filtered identities agree.')
    return counts
if __name__ == '__main__':
    small = check_profile((2, 1, 1), 38)
    check(small[9] == [0, 1470, 3756, 1470], 'small[9] == [0, 1470, 3756, 1470]')
    counts = check_profile((4, 2, 1), 54)
    r = 48
    coeffs = counts[r]
    polynomial = sp.Poly(sum((c * t ** j for j, c in enumerate(coeffs))), t)
    isolators = polynomial.intervals(eps=sp.Rational(1, 10 ** 16))
    roots = [float((a + b) / 2) for (a, b), m in isolators]
    check([int(sp.Rational(a + b, 2) * 10 ** 11) for (a, b), m in isolators] == [-500105643628, -200027360928, -100000000000, -49993160703, -19995775147, 0], '[int(sp.Rational(a + b, 2) * 10 ** 11) for (a, b), m in isolators] == [-500105643628, -200027360928, -100000000000, -499')
    print(f'r={r}, coefficients in ascending degree: {coeffs}')
    print('Figure 3, genus 20:', [f'{x:.14f}' for x in roots])
