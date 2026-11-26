# Replay receipt (v0.1.0-candidate, 26 November 2026)

Release replay by the publisher from fresh copies of the approved candidate code (results regenerated), working
directory = a copy of `code/`, macOS arm64. Interpreter as printed inside the run: CPython 3.12.11, sympy 1.14.0,
mpmath 1.3.0, numpy 2.5.5 (see [ENVIRONMENT.txt](ENVIRONMENT.txt)).

| Script | Command | Exit | Output vs recorded |
|---|---|---|---|
| verify.py | `uv run --no-project --python 3.12 --with sympy --with mpmath --with numpy python3 -B -O verify.py` | 0 | identical |
| check_third_residue.py | same | 0 | identical |
| verify_hurwitz.py | same | 0 | identical |
| verify_multiple.py | same | 0 | identical |
| check_hurwitz.py | same | 0 | identical |
| figure2.py | same | 0 | identical |
| occupancy.py | same | 0 | identical |
| census.py | same (about 50 s) | 0 | identical |
| Negative control | verify.py with the second-term factor `(D + profile.count(1) - 1)` on line 121 replaced by `(D + profile.count(1))`, under -O | 1 | `CHECK FAILED: sp.expand(Phi[D - 1] - second) == 0` |

"Identical" means byte for byte, including the final line `exit 0`. stderr was empty for every passing run. The lead's
earlier clean replay (24 November 2026, CPython 3.14.7) produced the recorded outputs. After the replay only the
title-page DOI line of manuscript.tex, the rebuilt PDF and the release records changed; code/, results/,
third-residue.tex and requirements.txt are byte-identical to the replayed candidate.
