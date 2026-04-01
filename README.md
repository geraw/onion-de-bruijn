# Onion De Bruijn

This repository is the code supplement for the paper:

Dor Genosar, Yotam Svoray, and Gera Weiss,
*Onion De Bruijn Sequences: Fixed-Window Counting by Growing the Alphabet*.

A PDF of the current paper draft is included at [`paper/Onion-De-Bruijn-Sequences.pdf`](paper/Onion-De-Bruijn-Sequences.pdf).

## What this repository provides

The repository contains reference implementations for the explicit arithmetic developed in the paper, together with a verification script that reproduces the finite computational checks discussed in the appendix.

- `order2_debruijn_arithmetic.py`
  Implements the order-2 rank and unrank maps, the canonical layer/offset coordinates, and direct addition and multiplication in the order-2 onion representation.
- `order3_debruijn_arithmetic.py`
  Implements the analogous order-3 rank and unrank maps, canonical coordinates, and direct addition and multiplication in the order-3 onion representation.
- `verify_results.py`
  Reproduces the finite computational checks from the paper, including De Bruijn/onion structure checks, counting formulas on small instances, and exhaustive arithmetic verification for the order-2 and order-3 constructions.

The code uses only the Python standard library.

## Requirements

- Python 3.10 or newer is recommended.

No external packages are required.

## Quick start

Clone the repository and run the verification script:

```bash
python verify_results.py
```

The script exits successfully if all checks pass and raises an exception if any tested identity fails.

## What is verified

The verification script checks the main explicit constructions on finite instances that are large enough to catch implementation mistakes while still being exhaustively tractable.

In particular, it verifies:

- the onion theorem on 14 small `(n, k)` instances
- the layer-count formula on 10 small instances
- the compatible onion-prefix count on 9 small instances
- universal structural and arithmetic consequences on 2032 exhaustively generated finite onion prefixes
- order-2 rank/unrank on 500 states and direct addition/multiplication on 14,400 pairs
- order-3 rank/unrank on 900 states and direct addition/multiplication on 8,100 pairs
- the exact carry-count formulas for both order-2 and order-3 direct arithmetic

These checks are intended as computational support for the explicit formulas and examples in the paper.

## Using the arithmetic code

The modules are written so they can be imported directly into a Python session or another script.

### Order 2

```python
from order2_debruijn_arithmetic import anti_rho2, rho2, add_states, mul_states

x = anti_rho2(10)
y = anti_rho2(7)

assert rho2(add_states(x, y)) == 17
assert rho2(mul_states(x, y)) == 70
```

### Order 3

```python
from order3_debruijn_arithmetic import anti_rho3, rho3, add_states3, mul_states3

x = anti_rho3(20)
y = anti_rho3(11)

assert rho3(add_states3(x, y)) == 31
assert rho3(mul_states3(x, y)) == 220
```

## Relation to the paper

The paper develops onion De Bruijn sequences as fixed-window counting systems and then gives explicit arithmetic for the reverse prefer-max onion sequences of orders 2 and 3.

This repository is meant to support three concrete uses:

- reproducing the finite computational checks reported in the paper
- experimenting with the order-2 and order-3 counting systems directly
- providing reference code for the rank/unrank maps, canonical coordinates, and direct arithmetic formulas

## Citation

If you use this code, please cite the accompanying paper:

```bibtex
@misc{genosar2026onion,
  author = {Dor Genosar and Yotam Svoray and Gera Weiss},
  title = {Onion De Bruijn Sequences: Fixed-Window Counting by Growing the Alphabet},
  year = {2026},
  note = {Manuscript},
  howpublished = {\url{https://github.com/geraw/onion-de-bruijn}}
}
```
