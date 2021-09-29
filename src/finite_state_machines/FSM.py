import itertools
from tqdm import tqdm


# Adapted from: https://stackoverflow.com/a/34325723/13526914
def printProgressBar(
    iteration, total, prefix="", suffix="", decimals=1, length=100, fill="█"
):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + "-" * (length - filledLength)
    print("\r%s |%s| %s%% %s" % (prefix, bar, percent, suffix), end="\r")


class FiniteStateMachine:
    """
    alphabet: set of length 1 strings
    num_states: integer
    start: integer
    accepting: set of integers
    transitions: dictionary, (integer, string): integer
            transition to -1 is explicit garbage state
    """

    def __init__(self, alphabet, num_states, start, accepting, transitions):
        self.alphabet = alphabet
        self.num_states = num_states
        self.start = start
        self.accepting = accepting
        self.transitions = transitions

        # It is assumed that there is no implicit garbage state if all possible
        # transitions are present.
        if self.num_states * len(self.alphabet) == len(self.transitions.keys()):
            self.explicit_garbage = True
        else:
            self.explicit_garbage = False

    def add_explicit_garbage(self):
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
                if (state, letter) not in self.transitions:
                    self.transitions[(state, letter)] = garbage_state
        for letter in self.alphabet:
            self.transitions[(garbage_state, letter)] = garbage_state
        self.explicit_garbage = True

    def brute_enumeration(self, length):
        """
        Returns the number of accepted words of a given length using brute force checking.
        You should probably never use this.
        """
        return len(self.words_generated(length))

    def words_generated(self, length):
        """
        Returns the list of accepted words of a given length using brute force checking.
        You should probably never use this.
        """
        possible_words = itertools.product(self.alphabet, repeat=length)
        return [w for w in possible_words if self.process_word(w)]

    def smart_enumeration(self, length, verbose=False):
        """
        Returns the counting sequence of accepted words up to a given length using
        dynamic programming.
        """
        enum = [1 if self.start in self.accepting else 0]

        counts = [0 for i in range(self.num_states)]
        counts[self.start] = 1

        lengths = (
            range(length)
            if not verbose
            else tqdm(range(length), dynamic_ncols=True, desc="Enumerating...")
        )
        for step in lengths:
            next_counts = [0 for i in range(self.num_states)]

            for (state, letter), next_state in self.transitions.items():
                next_counts[next_state] += counts[state]

            counts = list(next_counts)
            enum.append(
                sum([counts[i] for i in range(self.num_states) if i in self.accepting])
            )

        return enum

    def process_word(self, word):
        """
        Returns a boolean indicating whether the given word is accepted.
        """
        word = list(word)
        state = self.start
        for letter in word:
            if (state, letter) not in self.transitions:
                return False
            state = self.transitions[(state, letter)]
        return state in self.accepting

    def union(MA, MB, verbose=False):
        """
        Returns a new FiniteStateMachine object that is the union of MA and MB.
        The returned FSM has not been minimized.
        """
        return FiniteStateMachine.parallel(MA, MB, "union", verbose)

    def intersection(MA, MB, verbose=False):
        """
        Returns a new FiniteStateMachine object that is the intersection of MA and MB.
        The returned FSM has not been minimized.
        """
        return FiniteStateMachine.parallel(MA, MB, "intersection", verbose=verbose)

    def intersection_of_list(machines, verbose=False, minimize=True):
        """
        Computes the intersection of a given list with repeated pairwise
        intersections. It first intersects each consecutive pair, then the
        pairs of those, etc. This order seems faster in general.

        E.g.: (((A int B) int (C int D)) int ((E int F) int (G int H))) int ....
        """

        if verbose:
            print("\tIntersecting", len(machines), "FSAs")
        whichround = 1

        while len(machines) > 1:
            if verbose:
                print("\t\tRound", whichround, "with", len(machines), "machines")
            thisround = 1
            whichround += 1

            i = 0
            new_machines = []
            while len(machines) > 0:
                if len(machines) >= 2:
                    if verbose:
                        print("\t\t\t", thisround, "+", (thisround + 1))
                    thisround += 2

                    if minimize:
                        new_machines.append(
                            FiniteStateMachine.intersection(
                                machines[0], machines[1], verbose
                            ).minimize(inplace=False, verify=False, verbose=verbose)
                        )
                    else:
                        new_machines.append(
                            FiniteStateMachine.intersection(
                                machines[0], machines[1], verbose
                            )
                        )
                    machines = machines[2:]
                else:
                    if verbose:
                        print("\t\t\t", thisround, "alone")
                    thisround += 1
                    new_machines.append(machines[0])
                    machines = []
            machines = list(new_machines)

        return machines[0]

    def slower_intersection_of_list(machines, verbose=False, minimize=True):
        """
        Computes the intersection of a given list with by interesting in one
        new machine at a time. This seems to be the slower way to intersect
        many machines at once.

        E.g.: ((((A int B) int C) int D) int E) int ...
        """

        if verbose:
            print("\tIntersecting", len(machines), "FSAs")
        whichround = 1

        if len(machines) == 0:
            return -1
        elif len(machines) == 1:
            return machines[0]

        M = machines[0]
        for i in range(1, len(machines)):
            if verbose:
                print(
                    "\tIntersection with machine",
                    i,
                    "\t(",
                    M.num_states,
                    "+",
                    machines[i].num_states,
                    ")",
                )
            M = FiniteStateMachine.intersection(M, machines[i]).minimize(
                inplace=False, verify=True
            )

        return M

    def parallel(MA, MB, type, verbose=False):
        """
        This is mostly a helper function that returns a new machine that runs
        MA and MB in parallel. It is used to compute unions and intersections
        of machines.
        """
        MA.add_explicit_garbage()
        MB.add_explicit_garbage()
        alphabet = MA.alphabet.union(MB.alphabet)

        state_translation = []
        state_translation.append((MA.start, MB.start))
        state_translation_back = {(MA.start, MB.start): 0}

        start = 0

        states_to_process = set([start])
        states_processed = set()

        transitions = {}
        if verbose and (MA.num_states + MB.num_states > 100):
            printProgressBar(1, 1, prefix="buffer:", suffix="", length=50)
        max_so_far = 1

        while len(states_to_process) > 0:
            max_so_far = max(max_so_far, len(states_to_process))
            if verbose and (MA.num_states + MB.num_states > 100):
                printProgressBar(
                    len(states_to_process),
                    max_so_far,
                    prefix="buffer:",
                    suffix=str(len(states_to_process))
                    + "/"
                    + str(max_so_far)
                    + "  total: "
                    + str(len(state_translation))
                    + "      ",
                    length=50,
                )

            active_state = states_to_process.pop()
            active_state_meaning = state_translation[active_state]

            states_processed.add(active_state)

            for letter in alphabet:
                MA_dest = (
                    MA.transitions[(active_state_meaning[0], letter)]
                    if (active_state_meaning[0], letter) in MA.transitions
                    else -1
                )
                MB_dest = (
                    MB.transitions[(active_state_meaning[1], letter)]
                    if (active_state_meaning[1], letter) in MB.transitions
                    else -1
                )

                dest_state_meaning = (MA_dest, MB_dest)
                if (
                    dest_state_meaning not in state_translation
                    and dest_state_meaning != (-1, -1)
                ):
                    state_translation.append(dest_state_meaning)
                    state_translation_back[dest_state_meaning] = (
                        len(state_translation) - 1
                    )
                dest_state = state_translation_back[dest_state_meaning]

                # add to new transitions
                transitions[(active_state, letter)] = dest_state

                if (
                    dest_state not in states_to_process
                    and dest_state not in states_processed
                ):
                    states_to_process.add(dest_state)
        if verbose and (MA.num_states + MB.num_states > 100):
            print("\n\tResult has", len(state_translation), "states.")

        # now we need to set accepting states
        if type == "union":
            accepting = set(
                [
                    i
                    for i in range(len(state_translation))
                    if state_translation[i][0] in MA.accepting
                    or state_translation[i][1] in MB.accepting
                ]
            )
        elif type == "intersection":
            accepting = set(
                [
                    i
                    for i in range(len(state_translation))
                    if state_translation[i][0] in MA.accepting
                    and state_translation[i][1] in MB.accepting
                ]
            )
        else:
            raise Exception("invalid type")

        return FiniteStateMachine(
            alphabet, len(state_translation), start, accepting, transitions
        )

    @staticmethod
    def fsm_for_words_avoiding(word, alphabet=None):
        """
        Returns a FSM that accepts all words avoiding <word>. If alphabet=None,
        then the alphabet is inferred to be the set of symbols appearing in
        <word>.
        """

        # state 0 = no progress
        # state 1 = first letter matched
        #  ...
        # state len(word)+1 = matched

        if type(word) == str:
            word = list(word)
        if type(word) == tuple:
            word = list(word)

        if alphabet is None:
            alphabet = set(list(word))
        else:
            alphabet = set(alphabet)

        num_states = len(word) + 1
        start = 0
        accepting = set(range(len(word)))

        transitions = {}

        # if nothing is matched so far, begin trying to match
        for letter in alphabet:
            if letter == word[0]:
                transitions[(0, letter)] = 1
            else:
                transitions[(0, letter)] = 0

        # if we're in the middle of a match, we either move ahead or fall back
        for state in range(1, len(word)):
            word_seen = word[:state]
            for letter in alphabet:
                if letter == word[state]:
                    transitions[(state, letter)] = state + 1
                else:
                    # hard part, find the overlap
                    matching_against = word[:state] + [letter]
                    found_one = False
                    for new_start in range(1, len(matching_against)):
                        prefix = matching_against[new_start:]
                        if prefix == word[: len(prefix)]:
                            found_one = True
                            transitions[(state, letter)] = len(prefix)
                            break
                    if not found_one:
                        transitions[(state, letter)] = 0

        # if we're done with a match, stay in this state
        for letter in alphabet:
            transitions[(len(word), letter)] = len(word)

        return FiniteStateMachine(alphabet, num_states, start, accepting, transitions)

    # https://cs.stackexchange.com/questions/47061/help-with-understanding-hopcrofts-algorithm
    def minimize(self, inplace=True, verify=False, verbose=False):
        """
        Minimizes the FSM using Hopcroft's Algorithm.
            inplace: If True, modifies the calling object. If False, returns new object.
            verify: If True, performs a rough correctness check at the end.
            verbose: If True, prints intermediate information.
        """

        if not self.explicit_garbage:
            raise Exception("cannot minimize without explicit garbage")

        eq_classes = []
        if len(self.accepting) != 0:
            eq_classes.append(frozenset(self.accepting))
        if len(self.accepting) != self.num_states:
            eq_classes.append(
                frozenset(set(range(self.num_states)).difference(self.accepting))
            )
        eq_classes = set(eq_classes)

        processing = set([frozenset(self.accepting)])

        show_buffer = verbose

        if show_buffer:
            printProgressBar(1, 1, prefix="buffer:", suffix="", length=50)
        max_so_far = 1

        while len(processing) > 0:
            max_so_far = max(max_so_far, len(processing))
            if show_buffer:
                printProgressBar(
                    len(processing),
                    max_so_far,
                    prefix="     min buffer:",
                    suffix=str(len(processing))
                    + "/"
                    + str(max_so_far)
                    + "  total: "
                    + str(len(eq_classes))
                    + "      ",
                    length=50,
                )
            active_state = processing.pop()
            for active_letter in self.alphabet:
                states_that_move_into_active_state = frozenset(
                    [
                        state
                        for state in range(self.num_states)
                        if self.transitions[(state, active_letter)] in active_state
                    ]
                )

                copy_eq_classes = set(eq_classes)
                for checking_set in copy_eq_classes:
                    XintY = checking_set.intersection(
                        states_that_move_into_active_state
                    )
                    if len(XintY) == 0:
                        continue
                    XdiffY = checking_set.difference(states_that_move_into_active_state)
                    if len(XdiffY) == 0:
                        continue
                    eq_classes.remove(checking_set)
                    eq_classes.add(XintY)
                    eq_classes.add(XdiffY)
                    if checking_set in processing:
                        processing.remove(checking_set)
                        processing.add(XintY)
                        processing.add(XdiffY)
                    else:
                        if len(XintY) < len(XdiffY):
                            processing.add(XintY)
                        else:
                            processing.add(XdiffY)

        # now eq_classes are good to go, make them a list for ordering
        eq_classes = list(eq_classes)
        # print(eq_classes)
        if verbose:
            print()

        # need a backmap to prevent constant calls to index
        back_map = {}
        for i, eq in enumerate(eq_classes):
            for state in eq:
                back_map[state] = i

        new_alphabet = self.alphabet
        new_num_states = len(eq_classes)
        new_start = back_map[self.start]
        new_accepting = set([back_map[acc] for acc in self.accepting])

        new_transitions = {}
        for i, eq in enumerate(eq_classes):
            for letter in self.alphabet:
                new_transitions[(i, letter)] = back_map[
                    self.transitions[(list(eq)[0], letter)]
                ]

        if verify:
            T = FiniteStateMachine(
                new_alphabet, new_num_states, new_start, new_accepting, new_transitions
            )
            if T.smart_enumeration(20) != self.smart_enumeration(20):
                raise Exception("IncorrectMinimization")
            if verbose:
                print("\n\tminimization:", self.num_states, "->", new_num_states)

        if inplace:
            self.alphabet = new_alphabet
            self.num_states = new_num_states
            self.start = new_start
            self.accepting = new_accepting
            self.transitions = new_transitions
            return
        else:
            return FiniteStateMachine(
                new_alphabet, new_num_states, new_start, new_accepting, new_transitions
            )
