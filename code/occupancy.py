"""Exact check of the joint occupancy law (Proposition 7.2) in degrees 3, 4 and 5.

For every endpoint permutation, words are counted by dynamic programming in increasing upper letter,
recording the product, the orbit partition and the SET of upper letters used.
For each set S the transitive counts along one parity class are fitted exactly by an
exponential polynomial sum_k p_k(g) (k^2)^g with deg p_k < floor(d/(k+1)) (Lemma 2.3)
and checked on two further genera. The coefficient of D^r is compared with the
prediction  (2/d!) * [d in S] * prod_{b=2}^{d-1} w_b,  w_b = (b-1)/(D-b+1) if b in S else 1.
Standard library only. Failures raise SystemExit.
"""
from collections import Counter
from fractions import Fraction as Fr
from math import factorial


def check(condition, message=""):
    if not condition:
        raise SystemExit("CHECK FAILED: " + str(message))


def cycle_type(p):
    seen, out = set(), []
    for s in range(len(p)):
        if s not in seen:
            x, n = s, 0
            while x not in seen:
                seen.add(x)
                x = p[x]
                n += 1
            out.append(n)
    return tuple(sorted(out, reverse=True))


def words(d, limit):
    """counts[r][(perm, components, mask)] over monotone words of length r."""
    ident = tuple(range(d))
    allw = [Counter() for _ in range(limit + 1)]
    allw[0][ident, ident, 0] = 1
    for b in range(1, d):
        run = [Counter() for _ in range(limit + 1)]
        for length in range(1, limit + 1):
            for origin, new in ((allw, True), (run, False)):
                for (p, comp, mask), cnt in origin[length - 1].items():
                    for a in range(b):
                        q = list(p)
                        q[a], q[b] = q[b], q[a]
                        ca, cb = comp[a], comp[b]
                        joined = tuple(min(ca, cb) if x in (ca, cb) else x for x in comp)
                        run[length][tuple(q), joined, mask | (1 << b)] += cnt
        for length in range(1, limit + 1):
            allw[length].update(run[length])
    return allw


def solve(mat, rhs):
    n = len(mat)
    m = [row[:] + [v] for row, v in zip(mat, rhs)]
    for c in range(n):
        piv = next(i for i in range(c, n) if m[i][c] != 0)
        m[c], m[piv] = m[piv], m[c]
        inv = 1 / m[c][c]
        m[c] = [v * inv for v in m[c]]
        for i in range(n):
            if i != c and m[i][c] != 0:
                f = m[i][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[c])]
    return [row[-1] for row in m]


def main():
    total = 0
    for d, limit in ((3, 14), (4, 20), (5, 20)):
        D = d - 1
        allw = words(d, limit)
        cols = [(k, j) for k in range(1, d) for j in range(d // (k + 1))]
        n = len(cols)
        profiles = {}
        for r, states in enumerate(allw):
            for (p, comp, mask), cnt in states.items():
                if len(set(comp)) == 1:
                    profiles.setdefault(cycle_type(p), {}).setdefault(p, Counter())[r, mask] += cnt
        for mu, perms in sorted(profiles.items()):
            # occupancy counts are not central; the limit law is checked for EVERY endpoint
            r0 = d + len(mu) - 2
            lengths = list(range(r0, limit + 1, 2))
            check(len(lengths) >= n + 2, f"not enough lengths for d={d} {mu}")
            mat = [[Fr(g) ** j * Fr(k * k) ** g for k, j in cols] for g in range(n)]
            for perm, data in sorted(perms.items()):
                for mask in range(1 << d):
                    if mask & 1:
                        continue
                    vals = [Fr(data[r, mask]) for r in lengths]
                    sol = solve(mat, vals[:n])
                    for g in range(n + 2):
                        check(sum(c * Fr(g) ** j * Fr(k * k) ** g for c, (k, j) in zip(sol, cols)) == vals[g],
                              f"fit d={d} {perm} mask={mask}")
                    lead = sol[cols.index((D, 0))] / Fr(D) ** r0  # coefficient of D^r
                    pred = Fr(0)
                    if mask >> (d - 1) & 1:
                        pred = Fr(2, factorial(d))
                        for b in range(2, d):  # letter b has bit b-1
                            if mask >> (b - 1) & 1:
                                pred *= Fr(b - 1, D - b + 1)
                    check(lead == pred, f"occupancy d={d} {perm} mask={mask}: {lead} vs {pred}")
                    total += 1
        print(f"d={d}: joint occupancy leading coefficients agree for every endpoint permutation ({sum(len(v) for v in profiles.values())} permutations)")
    print(f"{total} endpoint/occupancy-set cases checked")


if __name__ == "__main__":
    main()
