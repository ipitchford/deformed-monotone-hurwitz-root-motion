"""Exact direct enumeration and finite-difference tests of the two pole amplitudes.
Uses only the Python standard library. Run with Python 3.10+.
"""
from collections import defaultdict
from fractions import Fraction
from math import factorial, prod

def check(condition, message=''):
    """Explicit check that also runs under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))

def partitions(n, maximum=None):
    if n == 0:
        yield ()
        return
    for first in range(min(n, maximum or n), 0, -1):
        for rest in partitions(n - first, first):
            yield ((first,) + rest)

def representative(profile):
    perm = []
    start = 0
    for size in profile:
        perm.extend([*range(start + 1, start + size), start])
        start += size
    return tuple(perm)

def next_state(state, a, b):
    perm, blocks = state
    new_perm = list(perm)
    new_perm[a], new_perm[b] = (new_perm[b], new_perm[a])
    lo, hi = sorted((blocks[a], blocks[b]))
    return (tuple(new_perm), tuple((lo if x == hi else x for x in blocks)))

def enumerate_connected(d, max_length, t):
    ident = tuple(range(d))
    table = [{(ident, ident): 1}] + [{} for _ in range(max_length)]
    for b in range(1, d):
        cache = {}
        active = {}
        updated = [table[0]]
        for length in range(1, max_length + 1):
            weights = dict(active)
            for state, weight in table[length - 1].items():
                weights[state] = weights.get(state, 0) + t * weight
            nxt = defaultdict(int)
            for state, weight in weights.items():
                if state not in cache:
                    cache[state] = [next_state(state, a, b) for a in range(b)]
                for new_state in cache[state]:
                    nxt[new_state] += weight
            active = nxt
            together = dict(active)
            for state, weight in table[length].items():
                together[state] = together.get(state, 0) + weight
            updated.append(together)
        table = updated
    connected = (0,) * d
    return {p: [row.get((representative(p), connected), 0) for row in table] for p in partitions(d)}

def residue(s, t):
    return t * prod((1 + Fraction(j * t, s - j) for j in range(1, s)))

def apply_shift_difference(seq, eigenvalue):
    return [b - eigenvalue * a for a, b in zip(seq, seq[1:])]

def main():
    scenarios = 0
    for d in (4, 5):
        m = d - 1
        max_length = 2 * (d * (m - 2) + 4) + 1
        for t in range(-2, d + 1):
            all_profiles = enumerate_connected(d, max_length, t)
            for profile, counts in all_profiles.items():
                lengths = list(range(2 + (d - len(profile)) % 2, max_length + 1, 2))
                f = profile.count(1)
                residual = [Fraction(factorial(d) * counts[L], 2) - m ** L * residue(m, t) + (m + f - 1) * (1 + t) * (m - 1) ** L * residue(m - 1, t) for L in lengths]
                for a in range(1, m - 1):
                    for _ in range(d):
                        residual = apply_shift_difference(residual, a * a)
                check(len(residual) >= 2 and (not any(residual)), (d, t, profile, residual))
                scenarios += 1
    print(f'Passed exact two-amplitude tests for {scenarios} degree/profile/marker combinations.')
    for t in (1, 2, 3):
        values = enumerate_connected(4, 9, t)[4,]
        print(f'd=4, profile=(4,), t={t}, counts:', values[3:10:2])
if __name__ == '__main__':
    main()
