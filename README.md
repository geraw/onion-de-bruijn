# Onion De Bruijn

This repository is the code supplement for the paper:

Dor Genosar, Yotam Svoray, and Gera Weiss,
*Onion De Bruijn Sequences: Fixed-Window Counting by Growing the Alphabet*.

The paper draft can be found in the arXiv [here](https://arxiv.org/abs/1906.06157).

## What this repository provides

The repository contains reference implementations for the explicit arithmetic developed in the paper, together with a verification script that reproduces the finite computational checks discussed in the appendix.

- `order2_debruijn_arithmetic.py`
  Implements the order-2 rank and unrank maps, the canonical layer/offset coordinates, direct addition and multiplication, and transported quotient/remainder in the order-2 onion representation.
- `order3_debruijn_arithmetic.py`
  Implements the analogous order-3 rank and unrank maps, canonical coordinates, direct addition and multiplication, and transported quotient/remainder in the order-3 onion representation.
- `tests/`
  Reproduces the finite computational checks from the paper, including De Bruijn/onion structure checks, counting formulas on small instances, and exhaustive arithmetic verification for the order-2 and order-3 constructions.
- `switching_activity.py`
  evaluates the bounded implementation experiment discussed in Section 4. It compares a binary counter modulo $9^4$, the reflected Gray encoding of that rank counter, and the order-$4$ onion counter truncated at maximal symbol $8$, generates the onion orbit by iterating the current successor rule from $0000$, checks that this orbit agrees exactly with the reverse prefer-max order on $[9]^4$, uses a moving-pointer realization of the onion state, and reports the resulting symbol-write, bit-toggle, and peak-to-average burst statistics for both binary and Gray-coded head pointers.

The code uses only the Python standard library and pytest for tests.

## Requirements

- Python 3.10 or newer is recommended.

No external packages are required.

## Quick start

Clone the repository, install pytest and run the tests:

```bash
python -m pytest tests
```

The tests will pass successfully if all checks pass and fail if any tested identity fails.

## What is verified

The tests check the main explicit constructions on finite instances that are large enough to catch implementation mistakes while still being exhaustively tractable.

In particular, it verifies:

- the onion theorem on 14 small `(n, k)` instances
- the layer-count formula on 10 small instances
- the compatible onion-prefix count on 9 small instances
- universal structural and arithmetic consequences on 2032 exhaustively generated finite onion prefixes
- order-2 rank/unrank on 1,296 states and direct addition and multiplication on 1,679,616 pairs
- order-3 rank/unrank on 4,096 states and direct addition and multiplication on 16,777,216 pairs
- the exact carry-count formulas for both order-2 and order-3 direct arithmetic
- the quotient and remainder layer bounds proved in the paper for transported division

These checks are intended as computational support for the explicit formulas and examples in the paper.

## Using the arithmetic code

The modules are written so they can be imported directly into a Python session or another script.

### Order 2

```python
from order2_debruijn_arithmetic import (
    add_states2,
    inverse_rho2,
    divmod_states2,
    mod_states2,
    mul_states2,
    rho2,
)

x = inverse_rho2(10)
y = inverse_rho2(7)

assert rho2(add_states2(x, y)) == 17
assert rho2(mul_states2(x, y)) == 70
assert divmod_states2(x, y) == (inverse_rho2(1), inverse_rho2(3))
assert rho2(mod_states2(x, y)) == 3
```

### Order 3

```python
from order3_debruijn_arithmetic import (
    add_states3,
    inverse_rho3,
    divmod_states3,
    mod_states3,
    mul_states3,
    rho3,
)

x = anti_rho3(20)
y = anti_rho3(11)

assert rho3(add_states3(x, y)) == 31
assert rho3(mul_states3(x, y)) == 220
assert divmod_states3(x, y) == (inverse_rho3(1), inverse_rho3(9))
assert rho3(mod_states3(x, y)) == 9
```

## Relation to the paper

The paper develops onion De Bruijn sequences as fixed-window counting systems and then gives explicit arithmetic for the reverse prefer-max onion sequences of orders 2 and 3, including direct addition and multiplication together with explicit quotient/remainder by rank transport.

This repository is meant to support three concrete uses:

- reproducing the finite computational checks reported in the paper
- experimenting with the order-2 and order-3 counting systems directly
- providing reference code for the rank/unrank maps, canonical coordinates, direct arithmetic formulas, and transported quotient/remainder

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
