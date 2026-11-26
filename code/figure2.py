"""Fixed-point identity (Proposition 4.5) and the comparison with CDM Figure 2 (Table 1).

1. Checks h_r(1,nu) = sum_i nu_i h_{r-1}(nu + e_i) (unnormalised form of
   H_g(1,nu) = sum_i (nu_i+1) H_g(nu+e_i)) against DIRECT word enumeration (verify.direct_words),
   for every nu of d-1 with d = 3..6 and every available length.
2. Computes the moving roots of every degree-7 profile at genus 20 by the CDM cut-and-join
   recursion and exact root isolation, checks them against the nine printed decimals of
   CDM Figure 2, and prints the first-order (two-term) and three-term predictions.
Requires SymPy and mpmath. Failures raise SystemExit.
"""
import os
import sys
from fractions import Fraction as Fr
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import sympy as sp
from verify import connected_single, cut_join, direct_words, integer_partitions


def check(condition, message=""):
    if not condition:
        raise SystemExit("CHECK FAILED: " + str(message))


def fixed_point_identity():
    cases = 0
    for d, limit in ((3, 14), (4, 14), (5, 13), (6, 11)):
        rows = connected_single(direct_words(d, limit), d)

        def h(profile, r):
            c = rows[r].get(tuple(sorted(profile, reverse=True)))
            return [c[j] for j in range(d)] if c else [0] * d

        for nu in integer_partitions(d - 1):
            for r in range(1, limit + 1):
                lhs = h((1,) + nu, r)
                rhs = [0] * d
                for i, part in enumerate(nu):
                    bigger = nu[:i] + (part + 1,) + nu[i + 1:]
                    for j, v in enumerate(h(bigger, r - 1)):
                        rhs[j] += part * v
                check(lhs == rhs, f"fixed-point identity d={d} nu={nu} r={r}")
                cases += 1
    print(f"fixed-point identity agrees with direct enumeration in {cases} (nu, length) cases, d=3..6")


CDM_FIG2 = {  # alpha_1, alpha_2, alpha_4, alpha_5 at genus 20, nine decimals, CDM Figure 2
    (7,): (-5.001267941, -2.000328365, -0.499917922, -0.199949295),
    (6, 1): (-5.001267941, -2.000328365, -0.499917922, -0.199949295),
    (5, 2): (-5.001056435, -2.000273606, -0.499931607, -0.199957751),
    (5, 1, 1): (-5.001232684, -2.000319238, -0.499920203, -0.199950704),
    (4, 3): (-5.001056438, -2.000273614, -0.499931605, -0.199957751),
    (4, 2, 1): (-5.001056436, -2.000273609, -0.499931607, -0.199957751),
    (4, 1, 1, 1): (-5.001173928, -2.000304027, -0.499924004, -0.199953053),
    (3, 3, 1): (-5.001056438, -2.000273614, -0.499931605, -0.199957751),
    (3, 2, 2): (-5.000880235, -2.000227985, -0.499943010, -0.199964796),
    (3, 2, 1, 1): (-5.001027065, -2.000266006, -0.499933507, -0.199958925),
    (3, 1, 1, 1, 1): (-5.001100492, -2.000285016, -0.499928756, -0.199955990),
    (2, 2, 2, 1): (-5.000880235, -2.000227985, -0.499943010, -0.199964796),
    (2, 2, 1, 1, 1): (-5.000978117, -2.000253332, -0.499936675, -0.199960882),
    (2, 1, 1, 1, 1, 1): (-5.001018906, -2.000263893, -0.499934035, -0.199959252),
    (1, 1, 1, 1, 1, 1, 1): (-5.001018906, -2.000263893, -0.499934035, -0.199959252),
}


def A(k, t):
    out = 1
    for j in range(1, k + 1):
        out *= t - 1 + mp.mpf(k) / j
    return out


def figure2():
    mp.mp.dps = 60
    d, g, D = 7, 20, 6
    T = sp.Symbol("t")
    worst1 = worst3 = 0
    print("profile | genus-20 roots (exact isolation) | rel. error of displacement: 2-term, 3-term")
    for profile, printed in CDM_FIG2.items():
        coeffs = cut_join(g, tuple(sorted(profile)))
        poly = sp.Poly(list(reversed(coeffs)), T)
        roots = sorted(mp.mpf(sp.Rational(a + b, 2)) for (a, b), m in
                       poly.intervals(eps=sp.Rational(1, 10 ** 40)))
        moving = [x for x in roots if x not in (0, -1) and abs(x + 1) > mp.mpf(10) ** -30 and abs(x) > mp.mpf(10) ** -30]
        check(len(moving) == 4, f"root count {profile}")
        for got, want in zip(moving, printed):
            check(abs(got - mp.mpf(want)) < mp.mpf("1.000001e-9"), f"Figure 2 (truncated to 9 decimals) mismatch {profile}: {got} vs {want}")
        f, c = profile.count(1), profile.count(2)
        r = 2 * g - 2 + d + len(profile)
        a = comb(d + f - 2, 2) - c
        Q = mp.mpf((d - 1) * (2 * c * (d - 2) + (d - 3) * (d - 4) + 2 * f * (d + f - 5))) / ((d - 4) * (d - 2))

        q1, q2 = mp.mpf(D - 1) / D, mp.mpf(D - 2) / D

        def two(t):  # scaled by D^-r
            return A(D, t) - (D + f - 1) * (1 + t) * q1 ** r * A(D - 1, t)

        def three(t):
            return two(t) + q2 ** r * A(D - 2, t) * (a * (1 + t) ** 2 - Q * t)

        errs = []
        for x, got in zip([1 - mp.mpf(D) / j for j in (1, 2, 4, 5)], moving):
            p2 = mp.findroot(two, got)
            p3 = mp.findroot(three, got)
            e2 = abs((p2 - got) / (got - x))
            e3 = abs((p3 - got) / (got - x))
            worst1, worst3 = max(worst1, e2), max(worst3, e3)
            errs.append((e2, e3))
        print(profile, [mp.nstr(x, 12) for x in moving],
              [f"{mp.nstr(e2, 2)}, {mp.nstr(e3, 2)}" for e2, e3 in errs])
    check(worst1 < mp.mpf("0.002") and worst3 < mp.mpf("1e-5"), "prediction accuracy")
    print(f"all 15 rows of CDM Figure 2 reproduced to the printed nine decimals; worst relative "
          f"displacement error: two-term {mp.nstr(worst1, 3)}, three-term {mp.nstr(worst3, 3)}")


if __name__ == "__main__":
    fixed_point_identity()
    figure2()
