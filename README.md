# FiniteStateMachines

This is a Python class to perform basic operations on finite state machines,
including union, intersection, and minimization.

Requires the `tqdm` library which can be inatlled via `pip`.

## Usage
```python
>>> from FSM import FiniteStateMachine as FSM

>>> M = FSM.fsm_for_words_avoiding("000", alphabet=["0","1"])
>>> M.smart_enumeration(10)
[1, 2, 4, 7, 13, 24, 44, 81, 149, 274, 504]

>>> N = FSM.fsm_for_words_avoiding("101", alphabet=["0","1"])
>>> N.smart_enumeration(10)
[1, 2, 4, 7, 12, 21, 37, 65, 114, 200, 351]

>>> M.intersection(N).smart_enumeration(10)
[1, 2, 4, 6, 9, 13, 19, 28, 41, 60, 88]

>>> M.union(N).smart_enumeration(10)
[1, 2, 4, 8, 16, 32, 62, 118, 222, 414, 767]
```

If this code was useful to you in your work, please consider citing it. A BibTeX entry is included below for your convenience.
```
@software{FiniteStateMachines,
  author       = {Jay Pantone},
  title        = {FiniteStateMachines: Version 1.0.0},
  month        = mar,
  year         = 2021,
  version      = {v1.0.0},
  doi          = {10.5281/zenodo.4592556},
  url          = {https://doi.org/10.5281/zenodo.4592556}
}
```

[![DOI](https://zenodo.org/badge/330728356.svg)](https://zenodo.org/badge/latestdoi/330728356)

Questions, comments, and improvements welcome!
