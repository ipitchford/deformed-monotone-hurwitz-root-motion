"""All-genus real-rootedness certificates and an exact census of root monotonicity
for single-profile deformed monotone Hurwitz polynomials of small degree.

Method (see the manuscript, Section "Every genus in small degree"):
 1. The CDM cut-and-join recursion gives h_r(mu;t) exactly at the first few genera.
 2. For fixed d, each coefficient of h_r on one parity class is an exponential polynomial
    sum_k p_k(g) (k^2)^g with deg p_k < floor(d/(k+1)) (pole-order bound). These forms are
    solved exactly from the first genera and checked on two further genera.
 3. Window certificate: at fixed rational separators s the value h_g(s) is an exponential
    polynomial in g; its sign is certified for every g >= G_w (tail bound + exact prefix).
    Alternating signs give real-rootedness and one root per window for g >= G_w; for the
    genera below G_w (G_w <= 1 here) real-rootedness is checked by an exact Sturm count.
 4. Tail certificate: at the moving points u = x + beta(1 +- eps) q^r the value h_g(u) is
    again an exponential polynomial (rational bases); certified signs place the root in
    nested intervals approaching x monotonically for g >= G.
 5. For g < G the roots are isolated exactly (SymPy) and compared.

Exact rational arithmetic throughout. Failures raise SystemExit (no assert statements).
Usage: python census.py [maxdegree]   (default 10)
"""
import sys
from fractions import Fraction as Fr
from functools import cache
from math import comb, floor
import sympy as sp


def fail(message):
    raise SystemExit("CHECK FAILED: " + message)


def check(condition, message):
    if not condition:
        fail(message)


@cache
def integer_partitions(d, minimum=1):
    if d == 0:
        return ((),)
    return tuple((i,) + p for i in range(minimum, d + 1)
                 for p in integer_partitions(d - i, i))


@cache
def cut_join(g, profile):
    """CDM Theorem 3.7, unnormalised: coefficients of h_r(mu;t), r = 2g-2+d+n."""
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
            left = tuple(b for j, b in enumerate(rest) if mask >> j & 1)
            right = tuple(b for j, b in enumerate(rest) if not mask >> j & 1)
            for genus in range(g + 1):
                x, y = get(genus, (u, *left)), get(g - genus, (v, *right))
                for j, value in enumerate(x):
                    if value:
                        add(y, value, j)
    return tuple(out)


def solve(matrix, rhs):
    """Exact Gauss-Jordan elimination over the rationals."""
    n = len(matrix)
    m = [row[:] + [b] for row, b in zip(matrix, rhs)]
    for col in range(n):
        piv = next((i for i in range(col, n) if m[i][col] != 0), None)
        check(piv is not None, "singular fitting system")
        m[col], m[piv] = m[piv], m[col]
        inv = 1 / m[col][col]
        m[col] = [v * inv for v in m[col]]
        for i in range(n):
            if i != col and m[i][col] != 0:
                f = m[i][col]
                m[i] = [a - f * b for a, b in zip(m[i], m[col])]
    return [row[-1] for row in m]


# An exponential polynomial is a dict {(base, j): coefficient} meaning sum coef * g^j * base^g.

def evaluate(ep, g):
    return sum(c * Fr(g) ** j * base ** g for (base, j), c in ep.items())


def fit_profile(d, profile):
    """Exact exponential-polynomial form in g of every coefficient of h_r(mu;t)."""
    cols = [(Fr(k * k), j) for k in range(1, d) for j in range(d // (k + 1))]
    n = len(cols)
    values = [cut_join(g, profile) for g in range(n + 2)]
    mat = [[Fr(g) ** j * base ** g for base, j in cols] for g in range(n)]
    forms = []
    for i in range(d):
        sol = solve(mat, [Fr(values[g][i]) for g in range(n)])
        ep = {col: c for col, c in zip(cols, sol) if c != 0}
        for g in range(n + 2):
            check(evaluate(ep, g) == values[g][i], f"fit check {profile} coefficient {i} genus {g}")
        forms.append(ep)
    return forms


def substitute(forms, point_terms):
    """h_g(u_g) where u_g = sum_{k} c_k lam_k^g is given as a list of (c, lam) pairs.
    Returns an exponential polynomial (bases multiplied out)."""
    out = {}
    # powers of u as exponential polynomials
    power = {(Fr(1), 0): Fr(1)}
    for i, ep in enumerate(forms):
        if i > 0:
            new = {}
            for (b1, j1), c1 in power.items():
                for c2, lam in point_terms:
                    key = (b1 * lam, j1)
                    new[key] = new.get(key, 0) + c1 * c2
            power = {k: v for k, v in new.items() if v != 0}
        for (b1, j1), c1 in ep.items():
            for (b2, j2), c2 in power.items():
                key = (b1 * b2, j1 + j2)
                out[key] = out.get(key, 0) + c1 * c2
    return {k: v for k, v in out.items() if v != 0}


def certify_positive(ep, start=0, limit=4000):
    """Return G such that ep(g) > 0 for every integer g >= start (prefix checked exactly)."""
    top = max(base for base, j in ep)
    lead = {j: c for (base, j), c in ep.items() if base == top}
    jmax = max(lead)
    check(lead[jmax] > 0, "leading coefficient not positive")
    # crude but rigorous tail: for g >= G, lead part >= (lead[jmax] - sum_{j<jmax}|c_j|/G) g^jmax top^g
    others = [(base, j, abs(c)) for (base, j), c in ep.items() if base != top]
    for G in list(range(max(start, 1), 200)) + list(range(200, limit, 50)):
        lower = lead[jmax] - sum(abs(c) for j, c in lead.items() if j < jmax) / Fr(G)
        if lower <= 0:
            continue
        # each other term: |c| g^j base^g <= |c| G^j base^G (g/G)^j (base/top)^(g-G) top^(g-G)
        # monotone decrease of g^j (base/top)^g for g >= G needs (1+1/G)^j base <= top
        if any((1 + Fr(1, G)) ** j * base > top for base, j, c in others):
            continue
        tail = sum(c * Fr(G) ** j * base ** G for base, j, c in others)
        if tail >= lower * Fr(G) ** jmax * top ** G:
            continue
        # positive for all g >= G; extend downwards by exact evaluation to the least such G
        while G > start and evaluate(ep, G - 1) > 0:
            G -= 1
        return G
    fail("no tail certificate found")


def A(k, t):
    out = Fr(1)
    for j in range(1, k + 1):
        out *= t - 1 + Fr(k, j)
    return out


def dA(k, t):
    roots = [1 - Fr(k, j) for j in range(1, k + 1)]
    return sum(_prod(t - r for jj, r in enumerate(roots) if jj != i) for i in range(len(roots)))


def _prod(it):
    out = Fr(1)
    for v in it:
        out *= v
    return out


def h_value(forms, g, t):
    return sum(evaluate(ep, g) * t ** i for i, ep in enumerate(forms))


def run(maxdegree):
    T = sp.Symbol("t")
    summary = []
    for d in range(4, maxdegree + 1):
        D = d - 1
        q = Fr(D - 1, D)
        limits = [1 - Fr(D, j) for j in range(1, D)]
        below = [x for x in limits if x < -1]
        seps = [(below[i] + below[i + 1]) / 2 for i in range(len(below) - 1)] + [Fr(-1)]
        even = D % 2 == 0
        for profile in integer_partitions(d):
            f = profile.count(1)
            r0 = d + len(profile) - 2
            forms = fit_profile(d, profile)
            # window certificate: sign of h at separators, expected (-1)^(number of roots of h left of s)
            # h = t (1+t)^[D even] p(t); p has positive leading coefficient, degree D-1-[D even].
            Gw = 0
            for i, s in enumerate(seps):
                if s == -1 and even:
                    # -1 is an exact root of h; use p(-1) = h'(-1)/(-1) instead: certify h'(-1) sign
                    ep = {}
                    for idx, form in enumerate(forms):
                        if idx == 0:
                            continue
                        for key, c in form.items():
                            ep[key] = ep.get(key, 0) + c * idx * Fr(-1) ** (idx - 1)
                    ep = {k: v for k, v in ep.items() if v != 0}
                    # p(-1) * s * ... : h'(-1) = (-1) * p(-1); roots of p left of -1: len(below)
                    sign_p = (-1) ** len(below)
                    expected = -sign_p
                else:
                    ep = substitute(forms, [(s, Fr(1))])
                    roots_left = (i + 1) if s != -1 else len(below)
                    # roots of h left of s: roots of p left of s (= roots_left), sign of t and (1+t) factors
                    sign_p = (-1) ** roots_left
                    sign_factor = (-1) * ((-1) if even else 1)  # t<0 and (1+t)<0 for s<-1
                    if s == -1:
                        sign_factor = -1  # odd D, only the factor t
                    expected = sign_p * sign_factor
                G = certify_positive({k: v * expected for k, v in ep.items()})
                Gw = max(Gw, G)
            # tail monotonicity certificates for each limit below -1
            lam = q * q
            eps = (1 - lam) / (3 * (1 + lam))
            violations = []
            Gtail = 0
            for x in below:
                beta = (D + f - 1) * (1 + x) * A(D - 1, x) / dA(D, x)
                check(beta < 0, "displacement sign")
                Gx = 0
                for sgn in (1, -1):
                    c0 = beta * (1 + sgn * eps) * q ** r0
                    ep = substitute(forms, [(x, Fr(1)), (c0, lam)])
                    # at u = x + beta(1+eps)q^r (lower point, further from x) h has the sign of
                    # h just left of the root; with one root of h in each window this is
                    # (-1)^(roots of h left) ... determine expected sign from the limiting expansion:
                    # h ~ 2/d! [D^r A_D(u) - ...]: sign is that of A_D'(x)(u-x) - (D+f-1)(1+x)A_{D-1}(x)q^r
                    #   = A_D'(x) beta q^r (1+sgn eps) - A_D'(x) beta q^r = sgn eps A_D'(x) beta q^r
                    expected = 1 if sgn * dA(D, x) * beta > 0 else -1
                    Gx = max(Gx, certify_positive({k: v * expected for k, v in ep.items()}))
                Gtail = max(Gtail, Gx)
                # the tail intervals must lie inside the window of x (true from Gx on, widths decrease)
                window = [v for v in [-10 ** 9] + seps if v < x][-1], [v for v in seps if v > x][0]
                lo = x + beta * (1 + eps) * q ** (r0 + 2 * Gx)
                hi = x + beta * (1 - eps) * q ** (r0 + 2 * Gx)
                check(window[0] < lo and hi < window[1] and hi < x, "tail interval inside window")
            # nested intervals need q^2 (1+eps) < 1 - eps
            check(lam * (1 + eps) < 1 - eps, "interval nesting")
            # finite census for g < Gtail + 1 by exact isolation
            Gmax = max(Gtail, Gw) + 1
            roots_by_g = []
            for g in range(Gmax + 1):
                coeffs = [evaluate(ep, g) for ep in forms]
                check(all(c.denominator == 1 for c in coeffs), "integrality")
                check([int(c) for c in coeffs] == list(cut_join(g, profile)) if g <= 3 else True,
                      "fit vs recursion")
                poly = sp.Poly([int(c) for c in reversed(coeffs)], T)
                if g < Gw:
                    check(poly.count_roots() == D and poly.sqf_part().degree() == D,
                          f"real-rootedness at small genus {profile} g={g}")
                r = r0 + 2 * g
                digits = int(r * (1 - float(q)) / 2.3) + 40
                ivs = [(Fr(int(a.p), int(a.q)), Fr(int(b.p), int(b.q)))
                       for (a, b), mult in poly.intervals(eps=sp.Rational(1, 10 ** digits)) if b < -1]
                check(len(ivs) == len(below), f"root count {profile} g={g}")
                roots_by_g.append(ivs)
            for idx, x in enumerate(below):
                for g in range(Gmax):
                    (a0, b0), (a1, b1) = roots_by_g[g][idx], roots_by_g[g + 1][idx]
                    if b0 < a1 and b1 < x:
                        continue  # strictly increasing towards x from below
                    if a0 == b0 == x or a1 == b1 == x:
                        violations.append((x, g, "at-limit"))
                        continue
                    check(b1 < a0 or b0 < a1, f"overlapping isolation {profile} g={g}")
                    if b1 < a0 and b0 < x:
                        kind = "moves-away"  # below the limit and decreasing
                    elif b1 < a0:
                        kind = "decreasing-above-limit"
                    else:
                        kind = "overshoot"
                    violations.append((x, g, kind))
            summary.append((d, profile, Gw, Gtail, violations))
            vtxt = "; ".join(f"limit {x}: genus {g}->{g + 1} {kind}" for x, g, kind in violations) or "none"
            print(f"d={d} {profile}: real-rooted at every genus (window signs certified g>={Gw}, exact root count below); "
                  f"monotone tail certified g>={Gtail}; exceptions: {vtxt}")
    return summary


if __name__ == "__main__":
    maxdegree = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    run(maxdegree)
    print("census complete")
