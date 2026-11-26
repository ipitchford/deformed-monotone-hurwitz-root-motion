# AI index: How the zeros of deformed monotone Hurwitz polynomials approach their limits (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex) with
[third-residue.tex](third-residue.tex); Section 1 is a plain introduction with the problem, an example, why it
matters, what was known, what is new and a roadmap).
DOI 10.5281/zenodo.24546686 (concept 10.5281/zenodo.24546685). Status: unrefereed candidate; internal AI-assisted
review and exact computer checks only.

## Normalisation

sigma in S_d of cycle type (profile) mu with n cycles; f = number of fixed points; D = d - 1. h_r(mu; t) = sum of
t^(hive number) over transitive monotone factorisations (a_1 b_1)...(a_r b_r) = sigma with a_i < b_i and
b_1 <= ... <= b_r; hive number = number of distinct b_i. Genus g from r = 2g - 2 + d + n. CDM = Coulter-Do-Moskovsky
(Canad. J. Math. 2025), whose deformed monotone Hurwitz numbers are these counts divided by the product of the parts.
A_k(t) = leading polynomial (eq. (A) in the paper), with zeros -(d-2), -(d-3)/2, ..., -1/(d-2), 0.

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A | h_r(mu;t) = (2/d!) [D^r A_D(t) - (D+f-1)(1+t) A_{D-1}(t) (D-1)^r + O((1+r)^d (D-2)^r)] for admissible r; the same for any finite list of profiles (marked version) | Theorem 3.1, Corollary 3.2 | written proof (Sections 2-3); leading term due to Accadia-Lewański-Ruzza; at t = 1 consistent with Yang; the fixed-point dependence parallels Yang's Theorem 1.3 for simple Hurwitz numbers |
| B | Each moving zero differs from its limit by an explicit multiple of ((D-1)/D)^r with the sign CDM predicted; eventual form of CDM Conjecture 3.15; large-genus cases of Conjectures 3.13-3.14 (credited to the leading term) | Theorem 4.1, Corollary 4.2 | written proof |
| C | Explicit threshold: r >= 200 d^2 log(20d) suffices (far from optimal; d = 7 needs r >= 48,429 while the census shows genus 0) | Theorem 4.4 | written proof |
| D | Removing a fixed point: identity from CDM's cut-and-join recursion explains the four exact coincidences in CDM's Figure 2; the two- and three-term expansions reproduce all 15 rows (worst relative displacement errors 1.7e-4 and 5.5e-10) | Proposition 4.5, Table 1 | written proof; figure2.py |
| E | Monotonicity cannot start at genus 0: (2,2) (degenerate: genus-0 zeros at the limits) and (2,2,2) (genuine reversal) | Theorem 5.1 | written proof with exact polynomials; verify.py |
| F | All 132 profiles of degrees 4-10: real-rooted at every genus; every zero moves strictly towards its limit at every step except (2,2), (2,2,2), (4,2,2), (4,4,2), (6,2,2), (4,2,2,2) (step 0->1) and (2,2,2,2), (2,2,2,2,2) (steps 0->1, 1->2); none in odd degree; never after genus 2 | Theorem 5.2 | exact certificate (census.py) with written justification; interlacing at small genus not certified |
| G | The coefficient of (D-j)^r depends only on the numbers of cycles of length <= j, for j < d/2 | Theorem 6.1 | written proof; verify.py (d = 6, 7, 8) |
| H | Third term explicit in fixed points and 2-cycles | Theorem 6.2 (third-residue.tex) | written proof; check_third_residue.py |
| I | Hive number is a sum of independent Bernoulli variables at large genus (mean (D+1)/2, variance (D^2-1)/(6D)); larger letter b used with probability (b-1)/(d-1), independently | Propositions 7.1-7.2 | written proof; occupancy.py |

## Evidence

| Check | Script | Output | Arithmetic |
|---|---|---|---|
| Both terms of claim A for every profile of degrees 2-5 up to length 24; pairs of profiles in degrees 3-4; cut-and-join recursion; Theorem 5.1 polynomials; all-genus certificates in degrees 4-5; short-cycle dependence (16/24/41 classes, d = 6/7/8) | [code/verify.py](code/verify.py) | [results/verify.out](results/verify.out) | exact (SymPy, Fractions) |
| Third term: symbolic identity in general d, f, c, t; 55 profiles of degrees 5-8 | [code/check_third_residue.py](code/check_third_residue.py) | [results/check_third_residue.out](results/check_third_residue.out) | exact |
| Independent standard-library checks: 91 profile/marker cases (two terms), 72 multi-profile cases | [code/verify_hurwitz.py](code/verify_hurwitz.py), [code/verify_multiple.py](code/verify_multiple.py) | [results/verify_hurwitz.out](results/verify_hurwitz.out), [results/verify_multiple.out](results/verify_multiple.out) | exact |
| Independent permutation-basis dynamic program; filtered identities for (2,1,1) and (4,2,1); genus-20 zeros of (4,2,1) (CDM Figure 3) | [code/check_hurwitz.py](code/check_hurwitz.py) | [results/check_hurwitz.out](results/check_hurwitz.out) | exact integers; roots isolated by SymPy, printed as interval midpoints |
| Proposition 4.5 (212 cases, degrees 3-6) and Table 1 | [code/figure2.py](code/figure2.py) | [results/figure2.out](results/figure2.out) | exact; mpmath 60 digits for root comparisons with CDM's decimals (tolerance 1e-9) |
| Proposition 7.2, every endpoint permutation in degrees 3-5 (2136 cases) | [code/occupancy.py](code/occupancy.py) | [results/occupancy.out](results/occupancy.out) | exact |
| Theorem 5.2 census, 132 profiles | [code/census.py](code/census.py) | [results/census.out](results/census.out) | exact exponential-polynomial sign certificates |

Numbers quoted on the release page: 132 profiles, degrees 4 to 10, 8 exceptional profiles, genus 2, probability
(b-1)/(d-1), limits -(d-2), ..., -1/(d-2), 4-cycle smallest zeros -2.618, -2.215, -2.086, -2.036 (genus 0-3), table
errors about 1.7e-4 and 5.5e-10, 15 rows of CDM's Figure 2.

## Limits

Fixed degree and growing genus only. The general threshold is far from optimal. Real zeros and the exact direction at
every genus are certified only up to d = 10; interlacing at small genus is not certified. The joint regime (degree and
genus both growing; cf. Barazer-Louf) and the b-deformed real-Grassmannian family (Coulter-Do Conjectures 4.21-4.22) are
open, as are the multi-profile direction where B(x_j) = 0 (Remark 4.3) and the rate of independence of the occupancy
indicators.

## Replay and records

Replay commands and the comparison policy (byte-for-byte equality with results/*.out, no tolerance):
[README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt), [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI
[.github/workflows/replay.yml](.github/workflows/replay.yml). Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art
and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography: [refs.bib](refs.bib). Licences:
[LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes: [MANIFEST.sha256](MANIFEST.sha256).
