import itertools

boolean_functions = {}
all_inputs = {}


def make_boolean_function(table):
    def boolean_function(input_set):
        return table[input_set]

    return boolean_function


def get_boolean_function(window_size, fi):
    global boolean_functions, all_outputs

    if window_size not in boolean_functions:
        boolean_functions[window_size] = {}

    if window_size not in all_inputs:
        all_inputs[window_size] = [
            x for x in itertools.product([0, 1], repeat=window_size)
        ]

    if fi not in boolean_functions[window_size]:
        all_outputs = itertools.product([0, 1], repeat=(2**window_size))
        for _ in range(fi):
            next(all_outputs)
        output_set = next(all_outputs)

        table = {}
        for idx, input_set in enumerate(all_inputs[window_size]):
            table[input_set] = output_set[idx]

        boolean_functions[window_size][fi] = make_boolean_function(table)

    return boolean_functions[window_size][fi]
