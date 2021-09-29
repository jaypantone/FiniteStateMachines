# FiniteStateMachines

This is a Python class to perform basic operations on finite state machines,
including union, intersection, and minimization.

It can be installed via pip with the command `pip install finite_state_machines`.

## Usage
```python
>>> from finite_state_machines import FiniteStateMachine as FSM

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

If this code was useful to you in your work, please consider citing it. To generate a BibTeX entry (or another format), click the badge below and locate the "Cite As" section.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4592556.svg)](https://zenodo.org/badge/latestdoi/330728356)


Questions, comments, and improvements welcome!
