"""
This module provides the WeightedFiniteStateMachine class.

HIGHLY EXPERIMENTAL! MANY FUNCTIONS HAVE NOT BEEN EXTENSIVELY TESTED.
"""
import itertools
from collections import defaultdict
from typing import Dict, List, Set, Tuple, Union

import sympy  # type: ignore

from .FSM import FiniteStateMachine


class WeightedFiniteStateMachine:  # pylint: disable=R0902
    """
    alphabet: set of length 1 strings
    num_states: integer
    start: integer
    accepting: set of integers
    weighted_transitions: dictionary, (integer, string): (integer, weight)
            where "weight" is a sympy polynomial in the variable "x"
            transition to -1 is explicit garbage state
    """

    def __init__(
        self,
        alphabet: Set[str],
        num_states: int,
        start: int,
        accepting: Set[int],
        weighted_transitions: Dict[
            Tuple[int, str], Tuple[int, Union[sympy.polys.polytools.Poly, sympy.Expr]]
        ],
    ):
        self.alphabet = alphabet
        self.num_states = num_states
        self.start = start
        self.accepting = accepting

        self.max_degree = 0
        self.x = sympy.symbols("x")
        self.weighted_transitions = dict()
        for ((state, letter), (next_state, weight)) in weighted_transitions.items():
            assert isinstance(
                weight, (sympy.Expr, sympy.polys.polytools.Poly)
            ), "weight must be a Sympy polynomial"
            weight_poly = weight.as_poly()
            assert weight_poly.free_symbols == {
                self.x
            }, "weights must be in the sympy variable 'x'"
            assert (
                weight_poly.coeff_monomial(1) == 0
            ), "weights must have a zero constant term"
            self.max_degree = max(self.max_degree, weight_poly.degree())
            self.weighted_transitions[(state, letter)] = (next_state, weight_poly)

        # keys are the size of words (by weight, not length!)
        # values are dicts whose
        #   keys represent which states they end in
        #   values are the sets of words of that size that end in that state
        self._word_cache: Dict[int, Dict[int, Set[str]]] = dict()

        # It is assumed that there is no implicit garbage state if all possible
        # transitions are present.
        if self.num_states * len(self.alphabet) == len(
            self.weighted_transitions.keys()
        ):
            self.explicit_garbage = True
        else:
            self.explicit_garbage = False

    def add_explicit_garbage(self) -> None:
        """
        Adds an explicit garbage state rather than relying on an implicit one.
        This is necessary for the minimize function to operate correctly.
        """
        if self.explicit_garbage:
            return
        garbage_state = self.num_states
        self.num_states += 1
        for state in range(self.num_states - 1):
            for letter in self.alphabet:
                if (state, letter) not in self.weighted_transitions:
                    self.weighted_transitions[(state, letter)] = (garbage_state, self.x)
        for letter in self.alphabet:
            self.weighted_transitions[(garbage_state, letter)] = (garbage_state, self.x)
        self.explicit_garbage = True

    def brute_enumeration(self, size: int) -> int:
        """
        Returns the number of accepted words of a given length using brute force
        checking. You should probably never use this.
        """
        return sum(1 for _ in self.words_generated(size))

    def words_generated(self, size: int) -> Set[str]:
        """
        Returns the set of words of a particular size accepted by the FSM.
        The result is cached in self.word_cache to speed up the generation to larger
        sizes.

        Only makes sense if all coefficients of all weight polynomials are 0 or 1, so
        raises an exception otherwise.
        """
        for (_, weight) in self.weighted_transitions.values():
            if not set(weight.all_coeffs()).issubset({0, 1}):
                raise NotImplementedError(
                    "words_generated only makes sense when all weight polynomials have "
                    "0/1 coefficients"
                )

        if size not in self._word_cache:
            if size == 0:
                self._word_cache[0] = {self.start: {""}}
            else:
                _ = self.words_generated(size - 1)

                new_words: Dict[int, Set[str]] = defaultdict(set)

                for prev_size in range(max(0, size - self.max_degree), size):
                    for (
                        (old_state, letter),
                        (new_state, weight),
                    ) in self.weighted_transitions.items():
                        if (
                            old_state in self._word_cache[prev_size]
                            and weight.coeff_monomial(self.x ** (size - prev_size)) == 1
                        ):
                            new_words[new_state].update(
                                w + letter
                                for w in self._word_cache[prev_size][old_state]
                            )
                self._word_cache[size] = new_words

        return set.union(
            *itertools.chain(
                self._word_cache[size][state]
                for state in self.accepting
                if state in self._word_cache[size]
            )
        )

    def smart_enumeration(self, size: int) -> List[int]:
        """
        Returns the counting sequence of accepted words up to a given size using
        dynamic programming.
        """
        enum = [1 if self.start in self.accepting else 0]

        counts = [
            [0 for i in range(self.num_states)] for j in range(self.max_degree + 1)
        ]
        counts[-1][self.start] = 1

        for _ in range(size):
            counts = counts[1:] + [[0 for i in range(self.num_states)]]

            for look_back in range(1, self.max_degree + 1):
                for (
                    (old_state, _),
                    (new_state, weight),
                ) in self.weighted_transitions.items():
                    counts[-1][new_state] += (
                        weight.coeff_monomial(self.x ** look_back)
                        * counts[-1 - look_back][old_state]
                    )

            enum.append(
                sum(
                    [
                        counts[-1][i]
                        for i in range(self.num_states)
                        if i in self.accepting
                    ]
                )
            )

        return enum

    def process_word(self, word_arg: str) -> bool:
        """
        Returns a boolean indicating whether the given word is accepted.
        Only works if all letters have weight 1.
        """
        print("WARNING: THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED.")
        if any(len(letter) != 1 for letter in self.alphabet):
            raise NotImplementedError(
                "Can only process words if all strings have length 1."
            )
        word = list(word_arg)
        state = self.start
        for letter in word:
            if (state, letter) not in self.weighted_transitions:
                return False
            (state, _) = self.weighted_transitions[(state, letter)]
        return state in self.accepting

    def union(
        self, other: "WeightedFiniteStateMachine", verbose: bool = False
    ) -> "WeightedFiniteStateMachine":
        """
        Returns a new WeightedFiniteStateMachine object that is the union of MA and MB.
        The returned WFSM has not been minimized.
        """
        print("WARNING: THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED.")
        fsms, back_dict = WeightedFiniteStateMachine.convert_WFSM_to_FSM([self, other])
        fsm_union = fsms[0].parallel(fsms[1], "union", verbose)
        return WeightedFiniteStateMachine.convert_FSM_to_WFSM(fsm_union, back_dict)

    def intersection(
        self, other: "WeightedFiniteStateMachine", verbose: bool = False
    ) -> "WeightedFiniteStateMachine":
        """
        Returns a new WeightedFiniteStateMachine object that is the intersection of MA
        and MB. The returned WFSM has not been minimized.
        """
        print("WARNING: THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED.")
        fsms, back_dict = WeightedFiniteStateMachine.convert_WFSM_to_FSM([self, other])
        fsm_int = fsms[0].parallel(fsms[1], "intersection", verbose)
        return WeightedFiniteStateMachine.convert_FSM_to_WFSM(fsm_int, back_dict)

    @staticmethod
    def intersection_of_list(
        machines: List["WeightedFiniteStateMachine"],
        verbose: bool = False,
        minimize: bool = True,
    ) -> "WeightedFiniteStateMachine":
        """
        Computes the intersection of a given list with repeated pairwise
        intersections. It first intersects each consecutive pair, then the
        pairs of those, etc. This order seems faster in general.

        E.g.: (((A int B) int (C int D)) int ((E int F) int (G int H))) int ....
        """
        print("WARNING: THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED.")
        fsms, back_dict = WeightedFiniteStateMachine.convert_WFSM_to_FSM(machines)
        fsm_int_of_list = FiniteStateMachine.intersection_of_list(
            fsms, verbose=verbose, minimize=minimize
        )
        return WeightedFiniteStateMachine.convert_FSM_to_WFSM(
            fsm_int_of_list, back_dict
        )

    # https://cs.stackexchange.com/questions/47061/help-with-understanding-hopcrofts-algorithm
    def minimize(
        self, verify: bool = False, verbose: bool = False
    ) -> "WeightedFiniteStateMachine":
        """
        Minimizes the FSM using Hopcroft's Algorithm.
            verify: If True, performs a rough correctness check at the end.
            verbose: If True, prints intermediate information.
        """
        print("WARNING: THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED.")
        fsms, back_dict = WeightedFiniteStateMachine.convert_WFSM_to_FSM([self])
        min_fsm = fsms[0].minimize(verify=verify, verbose=verbose)
        wfsm = WeightedFiniteStateMachine.convert_FSM_to_WFSM(min_fsm, back_dict)

        if verify:
            assert self.smart_enumeration(20) == wfsm.smart_enumeration(
                20
            ), "WFSM minimization failed verification"

        return wfsm

    @staticmethod
    def convert_WFSM_to_FSM(  # type: ignore
        machines: List["WeightedFiniteStateMachine"],
    ) -> Tuple[
        List[FiniteStateMachine], Dict[str, Tuple[str, sympy.polys.polytools.Poly]]
    ]:
        """
        In order to do operations like intersection, union, and minimize using the
        existing algorithms on the FSM object, we will convert the WFSMs to FSMs
        by extending the alphabet to assign each different (letter, weight)
        combination a new letter. When doing operations that involve more than one
        WFSM, we must do this in parallel so that the translations are the same.
        Returns a list of FSMs that can be operated on, as well as a dictionary to
        convert back, whose keys are "letters" like "TEMP-5" and whose values are the
        (letter, weight) pair that it should be reassigned to.
        """
        print("WARNING: THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED.")
        # build translation_dict and the reverse
        translation_dict: Dict[Tuple[str, sympy.polys.polytools.Poly], str] = dict()
        back_dict: Dict[str, Tuple[str, sympy.polys.polytools.Poly]] = dict()
        var_number = 0
        for machine in machines:
            for (_, letter), (_, weight) in machine.weighted_transitions.items():
                if (letter, weight) not in translation_dict:
                    var_name = f"TEMP-{var_number}"
                    translation_dict[(letter, weight)] = var_name
                    back_dict[var_name] = (letter, weight)
                    var_number += 1

        # build new FSMs
        # FiniteStateMachine(alphabet: Set[str], num_states: int, start: int,
        #           accepting: Set[int], transitions: Dict[Tuple[int, str], int])
        fsms = []
        for machine in machines:
            alphabet = set()
            num_states = machine.num_states
            start = machine.start
            accepting = machine.accepting
            transitions: Dict[Tuple[int, str], int] = dict()

            for (
                (state, letter),
                (new_state, weight),
            ) in machine.weighted_transitions.items():
                new_letter = translation_dict[(letter, weight)]
                alphabet.add(new_letter)
                transitions[(state, new_letter)] = new_state
            fsms.append(
                FiniteStateMachine(alphabet, num_states, start, accepting, transitions)
            )
        return fsms, back_dict

    @staticmethod
    def convert_FSM_to_WFSM(  # type: ignore
        machine: FiniteStateMachine,
        back_dict: Dict[str, Tuple[str, sympy.polys.polytools.Poly]],
    ) -> "WeightedFiniteStateMachine":
        """
        Converts one FSM back into a WFSM usng the back_dict produced by the
        original call to <convert_WFSM_to_FSM>.
        """
        print("WARNING: THIS FUNCTION HAS NOT BEEN EXTENSIVELY TESTED.")
        alphabet = set()
        num_states = machine.num_states
        start = machine.start
        accepting = machine.accepting
        weighted_transitions: Dict[
            Tuple[int, str], Tuple[int, Union[sympy.polys.polytools.Poly, sympy.Expr]]
        ] = dict()

        for (state, letter), new_state in machine.transitions.items():
            (real_letter, weight) = back_dict[letter]
            alphabet.add(real_letter)
            weighted_transitions[(state, real_letter)] = (new_state, weight)

        return WeightedFiniteStateMachine(
            alphabet, num_states, start, accepting, weighted_transitions
        )
