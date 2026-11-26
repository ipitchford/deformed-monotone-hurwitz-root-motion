"""Exact checks and all-genus certificate through degree five.

Run with Python 3.10+ and SymPy. No network or data files are used.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement, permutations
from math import factorial, prod
from functools import cache
import sympy as sp

def check(condition, message=''):
    """Explicit check that also runs under python -O."""
    if not condition:
        raise SystemExit('CHECK FAILED: ' + str(message))

def cycle_type(permutation):
    seen = set()
    lengths = []
    for start in range(len(permutation)):
        if start not in seen:
            x, size = (start, 0)
            while x not in seen:
                seen.add(x)
                x = permutation[x]
                size += 1
            lengths.append(size)
    return tuple(sorted(lengths, reverse=True))

def direct_words(d, limit):
    """Enumerate by DP over successive upper entries, preserving orbit partitions."""
    identity = tuple(range(d))
    all_words = [Counter() for _ in range(limit + 1)]
    all_words[0][identity, identity, 0] = 1
    for b in range(1, d):
        final_run = [Counter() for _ in range(limit + 1)]
        for length in range(1, limit + 1):
            for origin, new_upper in ((all_words, 1), (final_run, 0)):
                for (permutation, components, hive), count in origin[length - 1].items():
                    for a in range(b):
                        p = list(permutation)
                        p[a], p[b] = (p[b], p[a])
                        ca, cb = (components[a], components[b])
                        joined = tuple((min(ca, cb) if x in (ca, cb) else x for x in components))
                        key = (tuple(p), joined, hive + new_upper)
                        final_run[length][key] += count
        for length in range(1, limit + 1):
            all_words[length].update(final_run[length])
    return all_words

def connected_single(words, d):
    rows = []
    for states in words:
        grouped = defaultdict(Counter)
        for (permutation, orbit, hive), count in states.items():
            if len(set(orbit)) == 1:
                grouped[permutation][hive] += count
        by_profile = {}
        for permutation, counts in grouped.items():
            profile = cycle_type(permutation)
            check(by_profile.get(profile, counts) == counts, 'by_profile.get(profile, counts) == counts')
            by_profile[profile] = counts
        rows.append(by_profile)
    return rows

def fit(values, bases, var):
    """Fit p_k(g) k^(2g) given a rigorous pole-order bound."""
    cols = [(k, j) for k, degree in bases for j in range(degree)]
    m = len(cols)
    mat = sp.Matrix([[g ** j * k ** (2 * g) for k, j in cols] for g in range(m)])
    sol = mat.inv() * sp.Matrix(values[:m])
    result = defaultdict(lambda: 0)
    for c, (k, j) in zip(sol, cols):
        result[k * k] += c * var ** j
    for g, observed in enumerate(values):
        check(sum((sp.sympify(p).subs(var, g) * k ** g for k, p in result.items())) == observed, 'sum((sp.sympify(p).subs(var, g) * k ** g for k, p in result.items())) == observed')
    return result

def cross_difference(b, a, g):
    terms = defaultdict(lambda: 0)
    for i, bi in b.items():
        for j, aj in a.items():
            terms[i * j] += bi * aj.subs(g, g + 1) * j - bi.subs(g, g + 1) * aj * i
    return {k: sp.expand(p) for k, p in terms.items() if sp.expand(p) != 0}

def positive_on_nonnegative_integers(terms, g, start=0):
    """Certify strict positivity with a checked finite prefix and monotone tail bound."""
    largest = max(terms)
    c = terms[largest]
    check(c.is_Rational and c > 0, 'c.is_Rational and c > 0')
    others = {k: sp.Poly(p, g) for k, p in terms.items() if k < largest}
    for cut in range(1, 500):
        if not all(((cut + 1) ** p.degree() * k <= cut ** p.degree() * largest for k, p in others.items())):
            continue
        tail = sum((abs(coef) * cut ** j * k ** cut for k, p in others.items() for (j,), coef in p.terms()))
        if c * largest ** cut <= tail:
            continue
        check(all((sum((p.subs(g, i) * k ** i for k, p in terms.items())) > 0 for i in range(start, cut))), 'all((sum((p.subs(g, i) * k ** i for k, p in terms.items())) > 0 for i in range(start, cut)))')
        return cut
    raise ValueError('no certificate found')

def polynomial_A(k, t):
    return sp.expand(prod((t - 1 + sp.Rational(k, j) for j in range(1, k + 1))))

def verify_single(d, limit):
    rows = connected_single(direct_words(d, limit), d)
    profiles = sorted({p for row in rows for p in row})
    D = d - 1
    g, t = sp.symbols('g t')
    bases = [(k, d // (k + 1)) for k in range(1, d)]
    for profile in profiles:
        r0 = d + len(profile) - 2
        data = [rows[r].get(profile, Counter()) for r in range(r0, limit + 1, 2)]
        spectral = []
        for h in range(1, d):
            spectral.append(fit([s[h] for s in data], bases, g))
        Phi = {k: sp.expand(sum((s[k * k] * t ** i for i, s in enumerate(spectral, 1)))) for k in range(1, d)}
        top = sp.Rational(2, factorial(d)) * D ** r0 * polynomial_A(D, t)
        check(sp.expand(Phi[D] - top) == 0, 'sp.expand(Phi[D] - top) == 0')
        if d >= 3:
            second = -sp.Rational(2, factorial(d)) * (D - 1) ** r0 * (D + profile.count(1) - 1) * (1 + t) * polynomial_A(D - 1, t)
            check(sp.expand(Phi[D - 1] - second) == 0, 'sp.expand(Phi[D - 1] - second) == 0')
        for row in data:
            check(all((row[h] == row[d - h] for h in range(1, d))), 'all((row[h] == row[d - h] for h in range(1, d)))')
        if d in (4, 5):
            diff = cross_difference(spectral[1], spectral[0], g)
            start = 1 if profile == (2, 2) else 0
            cut = positive_on_nonnegative_integers(diff, g, start)
            if start:
                check(data[0] == Counter({1: 4, 2: 10, 3: 4}), 'data[0] == Counter({1: 4, 2: 10, 3: 4})')
                check(data[1] == Counter({1: 50, 2: 128, 3: 50}), 'data[1] == Counter({1: 50, 2: 128, 3: 50})')
            print(f'd={d}, profile={profile}: monotonicity from g={start}; tail bound g >= {cut}')
    print(f'd={d}: both polar coefficients checked for all {len(profiles)} profiles' if d > 2 else 'd=2: leading coefficient checked')
    return rows

def compose(p, q):
    return tuple((p[x] for x in q))

def inverse(p):
    return tuple((p.index(i) for i in range(len(p))))

def join_components(orbit, permutation):
    labels = list(orbit)
    for i, j in enumerate(permutation):
        a, b = (labels[i], labels[j])
        labels = [min(a, b) if v in (a, b) else v for v in labels]
    return tuple(labels)

def double_profile_check(d, limit=22):
    """Direct class average. The second branch permutation is determined by the first and word."""
    words = direct_words(d, limit)
    group = list(permutations(range(d)))
    sizes = Counter((cycle_type(p) for p in group))
    weights = {}
    rows = []
    for states in words:
        grouped = defaultdict(Counter)
        for (p, orbit, hive), count in states.items():
            key = (p, orbit)
            if key not in weights:
                w = Counter()
                for a in group:
                    b = compose(inverse(a), inverse(p))
                    if len(set(join_components(orbit, a))) == 1:
                        w[cycle_type(a), cycle_type(b)] += 1
                weights[key] = w
            for pair, multiplicity in weights[key].items():
                grouped[pair][hive] += count * multiplicity
        rows.append(grouped)
    D = d - 1
    g, t = sp.symbols('g t')
    bases = [(k, d // (k + 1)) for k in range(1, d)]
    pairs = list(combinations_with_replacement(sorted(sizes), 2))
    for a, b in pairs:
        r0 = 1 if (len(a) + len(b)) % 2 else 2
        pair, denom = ((a, b), sizes[a] * sizes[b])
        vals = [rows[r][pair] for r in range(r0, limit + 1, 2)]
        spectral = [fit([sp.Rational(v[h], denom) for v in vals], bases, g) for h in range(1, d)]

        def coeff(k):
            return sum((p[k * k] * t ** j for j, p in enumerate(spectral, 1)))
        L = (a.count(1) - 1) * (b.count(1) - 1)
        K = a.count(1) * b.count(1)
        B = 1 + L - K - (D + sp.Rational(L, D)) * t
        check(sp.expand(coeff(D) - sp.Rational(2, factorial(d)) * D ** r0 * polynomial_A(D, t)) == 0, 'sp.expand(coeff(D) - sp.Rational(2, factorial(d)) * D ** r0 * polynomial_A(D, t)) == 0')
        check(sp.expand(coeff(D - 1) - sp.Rational(2, factorial(d)) * (D - 1) ** r0 * polynomial_A(D - 1, t) * B) == 0, 'sp.expand(coeff(D - 1) - sp.Rational(2, factorial(d)) * (D - 1) ** r0 * polynomial_A(D - 1, t) * B) == 0')
    print(f'd={d}: both polar coefficients checked for all {len(pairs)} unordered double profiles')

@cache
def integer_partitions(d, minimum=1):
    if d == 0:
        return ((),)
    return tuple(((i,) + p for i in range(minimum, d + 1) for p in integer_partitions(d - i, i)))

@cache
def cut_join(g, profile):
    """CDM cut-and-join recursion, in the unnormalised h convention."""
    if g < 0:
        return ()
    if profile == (1,):
        return (int(g == 0),)
    a, *rest = profile
    out = [0] * sum(profile)

    def get(genus, values):
        return cut_join(genus, tuple(sorted(values)))

    def add(poly, weight=1, shift=0):
        for j, value in enumerate(poly):
            if value:
                out[j + shift] += weight * value
    for j, b in enumerate(rest):
        add(get(g, (a + b, *rest[:j], *rest[j + 1:])), b)
    if a > 1:
        short = get(g, (a - 1, *rest))
        add(short, -1)
        add(short, 1, 1)
    for u in range(1, a):
        v = a - u
        add(get(g - 1, (u, v, *rest)))
        for mask in range(1 << len(rest)):
            left = tuple((b for j, b in enumerate(rest) if mask >> j & 1))
            right = tuple((b for j, b in enumerate(rest) if not mask >> j & 1))
            for genus in range(g + 1):
                x, y = (get(genus, (u, *left)), get(g - genus, (v, *right)))
                for j, value in enumerate(x):
                    if value:
                        add(y, value, j)
    return tuple(out)

def cycle_count_check(d):
    """Compare all residues in the simple-pole range between profiles with equal short cycles."""
    g, t = sp.symbols('g t')
    bases = [(k, d // (k + 1)) for k in range(1, d)]
    columns = sum((order for k, order in bases))
    previous = {}
    comparisons = 0
    for profile in integer_partitions(d):
        vals = [cut_join(genus, profile) for genus in range(columns + 2)]
        r0 = d + len(profile) - 2
        spectrum = [fit([v[h] for v in vals], bases, g) for h in range(1, d)]
        for j in range((d - 1) // 2 + 1):
            pole = d - 1 - j
            coefficient = sp.expand(sp.Rational(factorial(d), 2 * pole ** r0) * sum((s[pole * pole] * t ** h for h, s in enumerate(spectrum, 1))))
            check(not coefficient.has(g), 'not coefficient.has(g)')
            quotient, remainder = sp.div(coefficient, polynomial_A(pole, t), t)
            check(remainder == 0 and sp.degree(quotient, t) <= j, 'remainder == 0 and sp.degree(quotient, t) <= j')
            key = (j, tuple((profile.count(a) for a in range(1, j + 1))))
            if key in previous:
                check(sp.expand(coefficient - previous[key]) == 0, 'sp.expand(coefficient - previous[key]) == 0')
                comparisons += 1
            previous[key] = coefficient
    print(f'd={d}: short-cycle dependence checked in {comparisons} repeated residue classes')
if __name__ == '__main__':
    for degree in range(2, 6):
        rows = verify_single(degree, 24)
        for r, profiles in enumerate(rows):
            for profile, counts in profiles.items():
                if r >= degree + len(profile) - 2:
                    g = (r - degree - len(profile) + 2) // 2
                    check(list(cut_join(g, tuple(sorted(profile)))) == [counts[j] for j in range(degree)], 'list(cut_join(g, tuple(sorted(profile)))) == [counts[j] for j in range(degree)]')
    rows = connected_single(direct_words(6, 9), 6)
    check([rows[7][2, 2, 2][j] for j in range(1, 6)] == [56, 424, 768, 424, 56], '[rows[7][2, 2, 2][j] for j in range(1, 6)] == [56, 424, 768, 424, 56]')
    check([rows[9][2, 2, 2][j] for j in range(1, 6)] == [2688, 19128, 34416, 19128, 2688], '[rows[9][2, 2, 2][j] for j in range(1, 6)] == [2688, 19128, 34416, 19128, 2688]')
    print('d=6: strict reversal checked by word counts at lengths 7 and 9')
    for degree in (3, 4):
        double_profile_check(degree)
    for degree in (6, 7, 8):
        cycle_count_check(degree)
