import itertools

boolean_functions = {}
all_inputs = {}


def make_boolean_function(table):
    def boolean_function(input_set):
        return table[input_set]

    return boolean_function


def nth_output_set(window_size, n):
    input_combinations = 2**window_size
    max_n = 2**input_combinations
    assert n < max_n
    bit = input_combinations
    output_set = []
    while bit > 0:
        if n % 2 == 1:
            output_set += [1]
        else:
            output_set += [0]
        n = n // 2
        bit -= 1
    return tuple(reversed(output_set))


def get_boolean_function(window_size, fi):
    global boolean_functions, all_outputs

    if window_size not in boolean_functions:
        boolean_functions[window_size] = {}

    if window_size not in all_inputs:
        all_inputs[window_size] = [
            x for x in itertools.product([0, 1], repeat=window_size)
        ]

    if fi not in boolean_functions[window_size]:
        output_set = nth_output_set(window_size, fi)

        table = {}
        for idx, input_set in enumerate(all_inputs[window_size]):
            table[input_set] = output_set[idx]

        boolean_functions[window_size][fi] = make_boolean_function(table)

    return boolean_functions[window_size][fi]
