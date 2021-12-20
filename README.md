# FiniteStateMachines

This package provides three classes for different types of finite state machines.
1. `FiniteStateMachine`, for classical finite state machines over an alphabet in which all letters are weighted equally;
2. `WeightedFiniteStateMachine`, for finite state machines over an alphabet in which each transition has a weight that is a polynomial in _x_;
3. `CombinatorialFSM`, for finite state machines with no alphabet at all, in which transitions between a pair of states are just recorded with a weight that is an expression in _x_ and possibly other variables.

It can be installed via pip with the command `pip install finite_state_machines`. Questions, comments, and improvements welcome!

---

## FiniteStateMachine

This is a Python class to perform basic operations on finite state machines, including union, intersection, and minimization.

```python
>>> from finite_state_machines import FiniteStateMachine as FSM

>>> M = FSM.fsm_for_words_avoiding("000", alphabet=["0","1"])
>>> M.enumeration(10)
[1, 2, 4, 7, 13, 24, 44, 81, 149, 274, 504]

>>> N = FSM.fsm_for_words_avoiding("101", alphabet=["0","1"])
>>> N.enumeration(10)
[1, 2, 4, 7, 12, 21, 37, 65, 114, 200, 351]

>>> M.intersection(N).words_generated(3)
{'001', '010', '011', '100', '110', '111'}

>>> M.intersection(N).enumeration(10)
[1, 2, 4, 6, 9, 13, 19, 28, 41, 60, 88]

>>> M.union(N).enumeration(10)
[1, 2, 4, 8, 16, 32, 62, 118, 222, 414, 767]
```

## WeightedFiniteStateMachine

In a weighted finite state machine, each transition is labeled with a weight that is a polynomial in _x_. This class has methods to generate words of a certain __size__ ("length" is no longer the correct metric) and to count words of a given size without generating them with a dynamic programming algorithm. Functions are provided to convert back and forth between weighted and non-weighted FSMs. Intersection and union of such machines are not precisely defined and thus not implemented, and there is currently no minimization function.

```python
>>> import sympy
>>> x = sympy.Symbol('x')
>>> from finite_state_machines import WeightedFiniteStateMachine as WFSM
>>> M = WFSM({"a", "b"}, 2, 0, {1}, {(0,"a"):(0,x),(0,"b"):(1,x),(1,"a"):(0,x**2),(1,"b"):(1,x**2)})
>>> [M.words_generated(i) for i in range(5)]
[set(), {'b'}, {'ab'}, {'aab', 'bb'}, {'aaab', 'abb', 'bab'}]

>>> M.enumeration(10)
[0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55]
```

## CombinatorialFSM

A combinatorial finite state machine (CFSM) has no alphabet. There are states (which in this implementation must be hashable), a single start state, a set of accepting states, and transitions between pairs of states that are weighted with sympy expressions that involve a main variable (_x_ be default) and possibly other variables.

Construction of a CFSM is a bit different from the previous two classes. After creating the class, you feed it transitions and weights. The weight is __added__ to any existing weight between these two same states.

The class comes with a minimization function and a function to write the transition matrix and solving routine to a Maple file.

```python
>>> import sympy
>>> x, y, C = sympy.symbols('x y C')
>>> from finite_state_machines import CombinatorialFSM
>>> M = CombinatorialFSM()
>>> M.add_transition(0, 1, x*y)
>>> M.add_transition(0, 1, x**2)
>>> M.add_transition(0, 2, x*y**2/C)
>>> M.add_transition(0, 3, x*y**2/C)
>>> M.add_transition(1, 0, x*(1+y/2))
>>> M.add_transition(2, 1, x**2)
>>> M.add_transition(3, 1, x**2)
>>> M.set_start(0)
>>> M.set_accepting({0, 1})
>>> M.enumeration(5)
[1,
 y,
 y**2/2 + y + 1,
 y**3/2 + y**2 + y/2 + 1 + 2*y**2/C,
 y**4/4 + y**3 + 2*y**2 + 2*y + y**3/C + 2*y**2/C,
 y**5/4 + y**4 + 3*y**3/2 + 2*y**2 + 5*y/2 + 1 + 2*y**4/C + 4*y**3/C]

>>> len(M.states)
4

>>> minimized = M.minimize()
>>> len(minimized.states)
3

>>> minimized.enumeration(5)
[1,
 y,
 y**2/2 + y + 1,
 y**3/2 + y**2 + y/2 + 1 + 2*y**2/C,
 y**4/4 + y**3 + 2*y**2 + 2*y + y**3/C + 2*y**2/C,
 y**5/4 + y**4 + 3*y**3/2 + 2*y**2 + 5*y/2 + 1 + 2*y**4/C + 4*y**3/C]

>>> with open('CFSM_example.txt', 'w') as f:
>>>     minimized.write_to_maple_file(f)
```

The produced Maple file is:
```
start := 1:
accepting := [1, 2]:
M := Matrix(3,3, storage=sparse):
M[1,2] := -x**2 - x*y:
M[1,3] := -2*x*y**2/C:
M[2,1] := x*(-y/2 - 1):
M[3,2] := -x**2:
M[1,1] := 1 + M[1,1]:
M[2,2] := 1 + M[2,2]:
M[3,3] := 1 + M[3,3]:
V := Vector(LinearAlgebra[Dimensions](M)[1]):
for a in accepting do V[a] := 1: od:
infolevel[solve] := 5:
xvec := LinearAlgebra[LinearSolve](M, V):
f := xvec[1]:
```

---

If this code is useful to you in your work, please consider citing it.

###### Bibtex entry:
```
@misc{FiniteStateMachines,
	author = {Jay Pantone},
	howpublished = {\url{https://github.com/jaypantone/FiniteStateMachines}},
	month = {September},
	note = {DOI: \url{https://doi.org/10.5281/zenodo.4592555}},
	title = {Finite{S}tate{M}achines},
	year = {2021}
}
```

###### Biblatex entry:
```
@software{FiniteStateMachines,
	author = {Jay Pantone},
	date = {2021},
	doi = {10.5281/zenodo.4592555},
	month = {9},
	title = {FiniteStateMachines},
	url = {https://github.com/jaypantone/FiniteStateMachines}
}
```
