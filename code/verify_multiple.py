"""Check the general residue against direct monotone walks and class tuples in S_4."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement, permutations, product
from math import factorial, prod

def check(condition, message=''):
    """Explicit check that also runs under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))

def multiply(p, q):
    return tuple((p[i] for i in q))

def join(blocks, a, b):
    low, high = sorted((blocks[a], blocks[b]))
    return tuple((low if i == high else i for i in blocks))

def profile(p):
    remaining = set(range(len(p)))
    lengths = []
    while remaining:
        j = next(iter(remaining))
        k = 0
        while j in remaining:
            remaining.remove(j)
            j = p[j]
            k += 1
        lengths.append(k)
    return tuple(sorted(lengths))

def endpoints(parts, classes, identity):
    modes = Counter()
    for tup in product(*(classes[p] for p in parts)):
        perm = blocks = identity
        for p in tup:
            perm = multiply(perm, p)
            for i, j in enumerate(p):
                blocks = join(blocks, i, j)
        inv = tuple((perm.index(i) for i in identity))
        modes[inv, blocks] += 1
    return modes

def word_layers(d, max_length, t):
    ident = tuple(range(d))
    current = {(ident, ident, -1): 1}
    tables = []
    for length in range(max_length + 1):
        sums = Counter()
        for (p, blocks, _), weight in current.items():
            sums[p, blocks] += weight
        tables.append(sums)
        nxt = defaultdict(int)
        for (p, blocks, prev), weight in current.items():
            for b in range(max(1, prev), d):
                for a in range(b):
                    h = list(p)
                    h[a], h[b] = (h[b], h[a])
                    nxt[tuple(h), join(blocks, a, b), b] += weight * (t if b != prev else 1)
        current = nxt
    return tables

def coefficients(parts, classes, tables, identity):
    modes = endpoints(parts, classes, identity)
    size = sum(modes.values())
    counts = []
    for layer in tables:
        total = 0
        for (p, orbits), multiplicity in modes.items():
            for (h, blocks), weight in layer.items():
                if p != h:
                    continue
                for i, rep in enumerate(orbits):
                    blocks = join(blocks, i, rep)
                if len(set(blocks)) == 1:
                    total += multiplicity * weight
        counts.append(Fraction(total, size))
    return counts

def R(s, t):
    return t * prod((1 + Fraction(j * t, s - j) for j in range(1, s)))

def main():
    d = 4
    ident = tuple(range(d))
    m = d - 1
    classes = defaultdict(list)
    for p in permutations(ident):
        classes[profile(p)].append(p)
    specifications = list(combinations_with_replacement(classes, 2))
    specifications += [((1, 1, 2), (1, 1, 2), (1, 3)), ((4,), (4,), (4,)), ((1, 1, 1, 1), (1, 1, 2), (1, 1, 2))]
    cases = 0
    for t in (-2, -1, 1, 2):
        tables = word_layers(d, 17, t)
        for parts in specifications:
            L0 = 2 + sum((d - len(p) for p in parts)) % 2
            n = len(parts)
            fixed = [p.count(1) for p in parts]
            A = Fraction(m) ** (2 - n) * prod((x - 1 for x in fixed))
            D = Fraction(d) ** (2 - n) * prod(fixed)
            B = 1 - m * t + A * (1 - Fraction(t, m)) - D
            counts = coefficients(parts, classes, tables, ident)
            residual = [Fraction(factorial(d), 2) * counts[L] - R(m, t) * m ** L - B * R(m - 1, t) * (m - 1) ** L for L in range(L0, len(tables), 2)]
            for _ in range(d):
                residual = [y - x for x, y in zip(residual, residual[1:])]
            check(residual and (not any(residual)), (parts, t, residual))
            cases += 1
    print(f'Passed {cases} profile-tuple/marker cases for the general residue.')
if __name__ == '__main__':
    main()
