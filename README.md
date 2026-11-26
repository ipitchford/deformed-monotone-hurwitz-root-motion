# How the zeros of deformed monotone Hurwitz polynomials approach their limits

Version 0.1.0-candidate, 26 November 2026. DOI [10.5281/zenodo.24546686](https://doi.org/10.5281/zenodo.24546686)
(all versions: [10.5281/zenodo.24546685](https://doi.org/10.5281/zenodo.24546685)). Author: Anonymous.
Release page: https://evidencepress.org/releases/deformed-monotone-hurwitz-root-motion/

**Status: unrefereed research candidate.** The results have written proofs that were checked by internal AI-assisted
review, with exact computer checks. No human referee has assessed them.

## What this is about

Take a shuffle (permutation) of d items and write it as a sequence of swaps of two items. Call the sequence
*monotone* if the larger item of each swap never decreases, and *transitive* if the swaps together connect all d
items. Count such sequences, giving each one a weight t for every distinct "larger item" it uses. The result is a
polynomial in t. Longer sequences correspond to higher *genus*.

Coulter, Do and Moskovsky (Canadian Journal of Mathematics, 2025) met these polynomials, which they call deformed
monotone Hurwitz numbers, in integrals over complex Grassmannians. Their data suggested that the zeros are all real,
that they interlace, and that, as the genus grows, each zero heads for a simple limit, -(d-2), -(d-3)/2, ..., -1/(d-2),
from a definite side (their Conjecture 3.15). For the 4-cycle, for example, the smallest zero is -2.618, -2.215, -2.086,
-2.036 at genus 0, 1, 2, 3, approaching -2 from below.

Why it matters: a polynomial with positive coefficients and only real zeros says that the number of larger items used
is exactly a sum of independent coin flips. Such hidden structure usually needs a mechanism.

What was known: the limits, and real zeros at very high genus, already follow from the leading large-genus term
found by Accadia, Lewański and Ruzza (2026). That term does not say from which side the zeros come.

## What is new

1. **The next term, exactly**, for every cycle type (and for several prescribed cycle types at once). It depends only
   on the number of fixed points, and its sign gives the direction Coulter, Do and Moskovsky predicted. So their
   Conjecture 3.15 holds in its eventual form, with an explicit, very generous threshold on the length.
2. **Not from the start.** For two 2-cycles and for three 2-cycles some zeros first move away from their limits.
3. **Every genus up to ten items.** For each of the 132 cycle types with 4 to 10 items, exact certificates prove real
   zeros at every genus and list exactly when a zero moves the wrong way: only 8 cycle types, all made of even cycles,
   and never after genus 2.
4. **Higher terms.** The third term depends only on fixed points and 2-cycles; early terms see only short cycles.
5. **Which larger items are used.** In a long random sequence, item b is used with probability (b-1)/(d-1),
   independently of the others.
6. **The published tables explained.** Some rows of the table of Coulter, Do and Moskovsky coincide exactly because a
   fixed point can be removed using their own recursion; the two- and three-term formulas reproduce all 15 rows of their
   Figure 2, with relative errors of about 1.7e-4 and 5.5e-10 in the displacement.

## Limits

All results are for fixed d and growing genus. The general threshold is far from optimal. Real zeros and the exact
direction at every genus are certified only up to d = 10, and interlacing at small genus is not certified. The regime
in which d and the genus grow together, and the b-deformed family from real Grassmannians, are open.

## Contents

| File | What it is |
|---|---|
| `manuscript.pdf`, `manuscript.tex`, `third-residue.tex` | The paper (17 pages); `third-residue.tex` is input by `manuscript.tex` |
| `code/verify.py` | Direct enumeration; both terms for every profile of degrees 2-5 up to length 24, pairs of profiles, cut-and-join recursion, Theorem 5.1 polynomials, short-cycle dependence |
| `code/check_third_residue.py` | The third term (symbolic identity and 55 profiles of degrees 5-8); imports `verify.py` |
| `code/verify_hurwitz.py`, `code/verify_multiple.py` | Independent standard-library checks (91 and 72 cases) |
| `code/check_hurwitz.py` | Independent dynamic program; genus-20 zeros of (4,2,1) |
| `code/figure2.py` | Fixed-point identity (212 cases) and Table 1 against CDM's Figure 2; imports `verify.py` |
| `code/occupancy.py` | Which larger items are used (2136 cases) |
| `code/census.py` | All 132 profiles of degrees 4-10 at every genus (about one minute) |
| `results/*.out` | Recorded outputs (each ends with the line `exit 0`); `results/python-version.txt` from the lead's replay |
| `refs.bib` | Bibliography |
| `AI_INDEX.md`, `RESEARCH_GATES.json`, `REVIEW_SUMMARY.md`, `REPLAY_RECEIPT.md`, `ENVIRONMENT.txt` | Records |

## Replay

Python 3.10 or later with SymPy, mpmath and NumPy (`requirements.txt`). From a directory holding all eight scripts:

    python3 -B -O verify.py        # likewise for the other seven scripts

Each script raises an error and exits 1 if any check fails, also under `python -O`. All arithmetic is exact
(Fractions, SymPy) except `figure2.py`'s 60-digit mpmath comparisons with published decimals. Comparison policy
(declared before tagging; also in CI, ENVIRONMENT.txt and REPLAY_RECEIPT.md): each output, followed by the line
`exit N`, equals `results/<script>.out` byte for byte; no tolerance is applied. GitHub CI
(`.github/workflows/replay.yml`) runs this on Linux, plus a negative control.

## Licence

Text, paper, outputs and records CC0-1.0; code MIT. See `LICENSE.md`. Cite with `CITATION.cff`.
