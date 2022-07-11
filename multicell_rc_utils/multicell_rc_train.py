import os
import sys
import struct
import math
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd
import seaborn as sns

from sklearn.linear_model import Lasso, LassoCV
from sklearn.model_selection import train_test_split
from sklearn.utils import parallel_backend

from vtk import vtkXMLPolyDataReader
from vtk.util.numpy_support import vtk_to_numpy

from .multicell_rc_params import Function
from .boolean_functions import get_boolean_function

STATES_FILE = "states"
BOOLEAN_FUNCTIONS_LASSO_MAX_ITERS = 10


def get_gene_values(output_file, num_genes, num_cells):
    """
    Extracts gene values for all cells from a simulation step file and returns them as a list.
    First num_genes values belong to cell 0, second num_genes to cell 1, etc.
    """
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(output_file)
    reader.Update()
    data = reader.GetOutput()

    assert data.GetNumberOfPoints() == num_cells

    cell_ids = [0] * num_cells
    for cell in range(num_cells):
        cell_ids[cell] = int(data.GetPointData().GetArray(1).GetTuple(cell)[0])

    gene_values = [0] * (num_cells * num_genes)
    for cell in range(num_cells):
        cell_id = cell_ids[cell]
        gene_count = 0
        for g in range(2, data.GetPointData().GetNumberOfArrays()):
            bits = data.GetPointData().GetArray(g).GetTuple(cell)[0]
            bits = struct.unpack("Q", struct.pack("d", bits))[0]
            width = min(num_genes - gene_count, 64)
            offset = cell_id * num_genes + gene_count
            vals = [int(n) for n in bin(bits)[2:].zfill(width)]
            gene_values[(offset) : (offset + width)] = vals[0:width]
            gene_count += width

    return gene_values


def get_states(
    num_genes, num_output_genes, num_cells, window_size, timesteps, output_dir, dump
):
    """
    Get gene values for all simulation steps as a list of lists.
    """

    states = [[] for _ in range(timesteps + 1)]
    for step in range(timesteps + 1):
        states[step] = get_gene_values(
            os.path.join(output_dir, "agent_0_0_0_" + str(step) + ".vtp"),
            num_genes,
            num_cells,
        )

    # Dump the states to a file
    if dump:
        with open(os.path.join(output_dir, STATES_FILE), "w") as f:
            for state in states:
                for i, bit in enumerate(state):
                    if i > 0:
                        f.write(" ")
                    f.write(str(bit))
                f.write("\n")

    return states


def get_cell_layers(num_cells, x_layers, y_layers, z_layers, cells_per_layer):
    """
    Returns a dictionary that maps cell id to cell layer.
    """
    cell_layer_map = {}
    skip_cells_per_layer = x_layers * y_layers - cells_per_layer
    cells_partial_layer_reach = math.ceil(num_cells / cells_per_layer)
    cells_full_layer_reach = math.floor(num_cells / cells_per_layer)
    cell_location = 0
    for cell in range(num_cells):
        if cell_location % (x_layers * y_layers) >= cells_per_layer:
            cell_location += skip_cells_per_layer
        elif (
            cells_partial_layer_reach < z_layers
            and cell_location // (x_layers * y_layers) == cells_full_layer_reach - 1
            and cell_location % (x_layers * y_layers) == cells_per_layer - 1
        ):
            cell_location += 1
        layer = cell_location // (x_layers * y_layers)
        cell_layer_map[cell] = layer
        cell_location += 1
    return cell_layer_map


def get_strains(output_file):
    """
    Returns a dictionary that maps cell id to strain.
    Takes Biocellion output as the input file.
    """
    strain_map = {}
    with open(output_file, "r") as f:
        for line in f:
            if line.startswith("Cell:"):
                tokens = line.split(",")
                cell = int(tokens[0].split(":")[1])
                strain = int(tokens[1].split(":")[1])
                strain_map[cell] = strain
    return strain_map


def process_output(
    input_signal_file,
    biocellion_output_file,
    output_dir,
    num_genes,
    num_cells,
    num_output_genes,
    num_output_cells,
    num_output_strains,
    output_cells_random,
    window_size,
    delay,
    timesteps,
    function,
    auxiliary_files,
    threads,
    warmup_steps,
    z_layers,
    y_layers,
    x_layers,
    input_signal_depth,
    cells_per_layer,
):
    """
    Builds the x and y lists for Lasso training. x consists of gene values from output cells and y is the ground truth for the given function and input signal. Also returns information about cell input signal reception.
    """
    with open(input_signal_file) as f:
        input_signal = [int(x) for x in f.readline().split()]

    states = get_states(
        num_genes,
        num_output_genes,
        num_cells,
        window_size,
        timesteps,
        output_dir,
        auxiliary_files,
    )
    states = states[
        1:
    ]  # We don't need initial state, only states 1+ are affected by the signal

    cell_layer_map = get_cell_layers(
        num_cells, x_layers, y_layers, z_layers, cells_per_layer
    )
    strain_map = get_strains(biocellion_output_file)

    cell_input_matches = [[] for _ in range(num_cells)]
    cell_input_constant = [[] for _ in range(num_cells)]
    for cell in range(num_cells):
        cell_input_matches[cell] = [False] * len(input_signal)
        cell_input_constant[cell] = [False] * len(input_signal)
        for i, (signal, genes) in enumerate(zip(input_signal, states)):
            cell_input_matches[cell][i] = signal == genes[cell * num_genes]
            cell_input_constant[cell][i] = (
                states[0][cell * num_genes] == genes[cell * num_genes]
            )
    input_signal_info = {}
    for layer in range(z_layers):
        input_signal_info[layer] = {
            "correct_cells": 0,
            "bad_cells": 0,
            "total_cells": 0,
        }

    for cell, (cell_matches, cell_constant) in enumerate(
        zip(cell_input_matches, cell_input_constant)
    ):
        layer = cell_layer_map[cell]
        if all(cell_matches):
            input_signal_info[layer]["correct_cells"] += 1
        elif not all(cell_constant):
            input_signal_info[layer]["bad_cells"] += 1
        input_signal_info[layer]["total_cells"] += 1
    cells_correct_input = 0
    cells_bad_input = 0
    for layer in range(z_layers):
        cells_correct_input += input_signal_info[layer]["correct_cells"]
        cells_bad_input += input_signal_info[layer]["bad_cells"]
    input_signal_info["correct_input"] = cells_correct_input
    input_signal_info["bad_input"] = cells_bad_input
    assert math.ceil(
        input_signal_info["correct_input"] / math.ceil(num_cells / z_layers)
    ) == min(input_signal_depth, z_layers)

    states = states[(window_size - 1) :]
    for state in states:
        assert len(state) == num_cells * num_genes

    if output_cells_random:
        output_cells = random.sample(range(num_cells), num_output_cells)
    else:
        output_cells = list(range(num_cells - 1, num_cells - num_output_cells - 1, -1))

    output_states = [[] for _ in range(len(states))]
    for i, state in enumerate(states):
        for cell in output_cells:
            strain = strain_map[cell]
            if strain >= num_output_strains:
                continue
            cell_state = state[(cell * num_genes) : ((cell + 1) * num_genes)]
            assert len(cell_state) == num_genes
            output_states[i] += cell_state[-num_output_genes:]
    for state in output_states:
        assert len(state) == len(output_states[0])
    not_output_strain = []
    for cell in output_cells:
        strain = strain_map[cell]
        if strain >= num_output_strains:
            not_output_strain.append(cell)
    for cell in not_output_strain:
        output_cells.remove(cell)
    if num_output_cells > 0 and len(output_cells) == 0:
        print(
            "\nNone of the output cells belong to any of the output strains. Exiting..."
        )
        sys.exit(1)
    x = output_states

    def calc_boolean_outputs(window_size, num_functions, input_signal, window_values):
        max_functions = 2 ** (2**window_size)
        assert num_functions <= max_functions
        fn_indices = random.sample(range(max_functions), num_functions)

        if num_functions == max_functions:
            fn_indices = sorted(fn_indices)

        outputs = [
            ([0] * (len(input_signal) + 1 - window_size)) for _ in range(num_functions)
        ]
        for fi in range(num_functions):
            f = get_boolean_function(window_size, fn_indices[fi])
            f_output = 0
            for i, j in enumerate(range(window_size, len(input_signal) + 1)):
                bits = window_values(f_output, input_signal, j)
                f_output = f(bits)
                outputs[fi][i] = f_output
            assert len(input_signal) == len(outputs[fi]) + window_size - 1
        return outputs, fn_indices

    fn_indices = None
    if function == Function.PARITY.value:
        parity = [0] * (len(input_signal) + 1 - window_size)
        for i, j in enumerate(range(window_size, len(input_signal) + 1)):
            bitsum = 0
            for bit in input_signal[(j - window_size) : j]:
                bitsum += bit
            if bitsum % 2 == 0:
                parity[i] = 0
            else:
                parity[i] = 1
        assert len(input_signal) == len(parity) + window_size - 1
        y = parity
    elif function == Function.MEDIAN.value:
        median = [0] * (len(input_signal) + 1 - window_size)
        for i, j in enumerate(range(window_size, len(input_signal) + 1)):
            bitsum = 0
            for bit in input_signal[(j - window_size) : j]:
                bitsum += bit
            if bitsum > window_size // 2:
                median[i] = 1
            else:
                median[i] = 0
        assert len(input_signal) == len(median) + window_size - 1
        y = median
    elif function == Function.THREE_BIT.value or function == Function.FIVE_BIT.value:
        y, fn_indices = calc_boolean_outputs(
            window_size,
            256,
            input_signal,
            lambda f_output, input_signal, j: tuple(
                input_signal[(j - window_size) : j]
            ),
        )
    elif (
        function == Function.THREE_BIT_RECURSIVE.value
        or function == Function.FIVE_BIT_RECURSIVE.value
    ):
        y, fn_indices = calc_boolean_outputs(
            window_size,
            256,
            input_signal,
            lambda f_output, input_signal, j: (f_output,)
            + tuple(input_signal[(j - window_size + 1) : j]),
        )
    else:
        assert False, f"Invalid function: {function}"

    if delay > 0:
        x = x[delay:]
        # For recursive functions, we get a list of arrays instead of a single array
        if type(y[0]) == list:
            for i in range(len(y)):
                y[i] = y[i][:-delay]
        else:
            y = y[:-delay]

    if warmup_steps > 0:
        x = x[warmup_steps:]
        if type(y[0]) == list:
            for i in range(len(y)):
                y[i] = y[i][warmup_steps:]
        else:
            y = y[warmup_steps:]

    if type(y[0]) == list:
        for i in range(len(y)):
            assert len(x) == len(y[i])
    else:
        assert len(x) == len(y)

    for s in x:
        assert len(s) == len(output_cells) * num_output_genes

    return x, y, input_signal_info, output_cells, fn_indices


def run_lasso(x, y, threads, max_iter=1000):
    x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.25)

    lasso = LassoCV(n_jobs=threads, selection="random", tol=0.05, max_iter=max_iter)
    lasso.fit(x_train, y_train)

    train_predicted = [1 if x > 0.5 else 0 for x in lasso.predict(x_train)]
    test_predicted = [1 if x > 0.5 else 0 for x in lasso.predict(x_test)]

    train_score = lasso.score(x_train, y_train)
    test_score = lasso.score(x_test, y_test)

    train_accuracy = sum(
        [1 if a == b else 0 for a, b in zip(train_predicted, y_train)]
    ) / len(train_predicted)
    test_accuracy = sum(
        [1 if a == b else 0 for a, b in zip(test_predicted, y_test)]
    ) / len(test_predicted)

    return train_accuracy, test_accuracy, lasso


def train_lasso(
    x,
    y,
    num_cells,
    num_output_cells,
    num_output_strains,
    num_output_genes,
    biocellion_output_file,
    output_dir,
    threads,
    output_cells,
    x_layers,
    y_layers,
    z_layers,
    cells_per_layer,
):
    """
    Train Lasso on the provided x and y lists and plots a histogram of features used from each cell and strain.
    """

    boolean_functions_train_accuracies = None
    boolean_functions_test_accuracies = None
    boolean_functions_lassos = None

    # In case of recursive functions, we have a list of functions instead of a single function
    if type(y[0]) == list:
        boolean_functions_train_accuracies = []
        boolean_functions_test_accuracies = []
        boolean_functions_lassos = []
        for i in range(len(y)):
            train_accuracy, test_accuracy, lasso = run_lasso(
                x, y[i], threads, BOOLEAN_FUNCTIONS_LASSO_MAX_ITERS
            )
            boolean_functions_train_accuracies.append(train_accuracy)
            boolean_functions_test_accuracies.append(test_accuracy)
            boolean_functions_lassos.append(lasso)
            if (i % 10 == 0) and (i != 0):
                print(f"LASSO finished for {i}/{len(y)}")
        train_accuracy = sum(boolean_functions_train_accuracies) / len(
            boolean_functions_train_accuracies
        )
        test_accuracy = sum(boolean_functions_test_accuracies) / len(
            boolean_functions_test_accuracies
        )
    else:
        train_accuracy, test_accuracy, lasso = run_lasso(x, y, threads)

    cell_layer_map = get_cell_layers(
        num_cells, x_layers, y_layers, z_layers, cells_per_layer
    )
    strain_map = get_strains(biocellion_output_file)

    output_cells_rank = {}
    for rank, cell in enumerate(output_cells):
        output_cells_rank[cell] = rank

    if type(y[0]) == list:
        cell_feature_data = pd.DataFrame()
        cell_weights = None
        strain_feature_data = pd.DataFrame()
    else:
        cell_feature_data = {"index": [], "cell": [], "features": [], "layer": []}
        cell_weights = {}
        strain_feature_data = {"index": [], "strain": [], "features": []}
        for cell in output_cells:
            cell_layer = cell_layer_map[cell]
            strain = strain_map[cell]
            cell_rank = output_cells_rank[cell]
            if strain >= num_output_strains:
                assert False
            coeff_start = cell_rank * num_output_genes
            coeff_end = (cell_rank + 1) * num_output_genes
            cell_weights[cell] = []
            for weight in lasso.coef_[coeff_start:coeff_end]:
                cell_weights[cell].append(weight)
            feature_count = np.sum(lasso.coef_[coeff_start:coeff_end] != 0)
            if feature_count > 0:
                cell_feature_data["cell"].append(cell)
                cell_feature_data["features"].append(feature_count)
                cell_feature_data["layer"].append(cell_layer)
                if strain in strain_feature_data["strain"]:
                    idx = strain_feature_data["strain"].index(strain)
                    strain_feature_data["features"][idx] += feature_count
                else:
                    strain_feature_data["strain"].append(strain)
                    strain_feature_data["features"].append(feature_count)
        cell_feature_data["index"] = [
            num for num in range(len(cell_feature_data["cell"]))
        ]
        strain_feature_data["index"] = [
            num for num in range(len(strain_feature_data["strain"]))
        ]

        cell_feature_data = pd.DataFrame(cell_feature_data)
        strain_feature_data = pd.DataFrame(strain_feature_data)

        if len(cell_feature_data) > 0:
            plt.figure()
            sns.barplot(
                data=cell_feature_data,
                x="index",
                y="features",
                hue="layer",
                dodge=False,
            )
            plt.xticks([])
            plt.xlabel("Cell")
            plt.ylabel("Features", rotation=0, labelpad=30)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "cell_feature_count.png"))

        if len(strain_feature_data) > 0:
            plt.figure()
            sns.barplot(data=strain_feature_data, x="index", y="features", dodge=False)
            plt.xticks([])
            plt.xlabel("Strain")
            plt.ylabel("Features", rotation=0, labelpad=30)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, "strain_feature_count.png"))

    max_features = len(x[0])
    return (
        train_accuracy,
        test_accuracy,
        cell_feature_data,
        strain_feature_data,
        max_features,
        cell_weights,
        boolean_functions_train_accuracies,
        boolean_functions_test_accuracies,
    )
