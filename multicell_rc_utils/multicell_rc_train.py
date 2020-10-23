import os
import sys
import shutil
import subprocess
import struct
import math
import matplotlib.pyplot as plt
import numpy as np
import random

from vtk import vtkXMLPolyDataReader
from vtk.util.numpy_support import vtk_to_numpy

STATES_FILE = 'states'

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

  cell_ids = [ 0 ] * num_cells
  for cell in range(num_cells):
    cell_ids[cell] = int(data.GetPointData().GetArray(1).GetTuple(cell)[0])

  gene_values = [ 0 ] * (num_cells * num_genes)
  for cell in range(num_cells):
    cell_id = cell_ids[cell]
    gene_count = 0
    for g in range(2, data.GetPointData().GetNumberOfArrays()):
      bits = data.GetPointData().GetArray(g).GetTuple(cell)[0]
      bits = struct.unpack('Q', struct.pack('d', bits))[0]
      width = min(num_genes - gene_count, 64)
      offset = cell_id * num_genes + gene_count
      vals = [ int(n) for n in bin(bits)[2:].zfill(width) ]
      gene_values[(offset):(offset + width)] = vals[0:width]
      gene_count += width

  return gene_values

def get_states(num_genes, num_output_genes, num_cells, window_size, timesteps, output_dir, dump):
  """
    Get gene values for all simulation steps as a list of lists.
  """

  states = [ [] for _ in range(timesteps + 1) ]
  for step in range(timesteps + 1):
    states[step] = get_gene_values(os.path.join(output_dir, 'agent_0_0_0_' + str(step) + '.vtp'), num_genes, num_cells)

  # Dump the states to a file
  if dump:
    with open(os.path.join(output_dir, STATES_FILE), 'w') as f:
      for state in states:
        for i, bit in enumerate(state):
          if i > 0:
            f.write(' ')
          f.write(str(bit))
        f.write('\n')

  return states

def get_cell_types(output_file):
  """
    Returns a dictionary that maps cell id to cell type.
    Takes Biocellion output as the input file.
  """
  cell_type_map = {}
  with open(output_file, 'r') as f:
    for line in f:
      if line.startswith("Cell:"):
        tokens = line.split(',')
        cell = int(tokens[0].split(':')[1])
        cell_type = int(tokens[1].split(':')[1])
        cell_type_map[cell] = cell_type
  return cell_type_map

def process_output(input_signal_file, biocellion_output_file, output_dir, num_genes, num_cells, num_output_genes, num_output_cells, num_output_cell_types, output_cells_random, window_size, delay, timesteps, function, auxiliary_files, threads, warmup_steps, z_layers, y_layers, x_layers, input_signal_depth):
  """
    Builds the x and y lists for Lasso training. x consists of gene values from output cells and
    y is the ground truth for the given function and input signal.
    Also returns information about cell input signal reception.
  """
  with open(input_signal_file) as f:
    input_signal = [ int(x) for x in f.readline().split() ]

  states = get_states(num_genes, num_output_genes, num_cells, window_size, timesteps, output_dir, auxiliary_files)

  cell_input_matches = [ [] for _ in range(num_cells) ]
  cell_input_constant = [ [] for _ in range(num_cells) ]
  for cell in range(num_cells):
    cell_input_matches[cell] = [ False ] * len(input_signal)
    cell_input_constant[cell] = [ False ] * len(input_signal)
    for i, (signal, genes) in enumerate(zip(input_signal, states)):
      cell_input_matches[cell][i] = (signal == genes[cell * num_genes])
      cell_input_constant[cell][i] = (states[0][cell * num_genes] == genes[cell * num_genes])
  input_signal_info = {}
  for layer in range(z_layers):
    input_signal_info[layer] = { "correct_cells" : 0, "bad_cells" : 0, "total_cells" : 0 }

  for cell, (cell_matches, cell_constant) in enumerate(zip(cell_input_matches, cell_input_constant)):
    layer = cell // (x_layers * y_layers)
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
  assert round(input_signal_info["correct_input"] / (num_cells / z_layers)) == input_signal_depth

  states = states[(window_size - 1):]
  for state in states:
    assert len(state) == num_cells * num_genes

  cell_types_map = get_cell_types(biocellion_output_file)

  if output_cells_random:
    output_cells = random.sample(range(num_cells), num_output_cells)
  else:
    output_cells = range(num_cells - 1, num_cells - num_output_cells - 1, -1)

  output_states = [ [] for _ in range(len(states)) ]
  for i, state in enumerate(states):
    for cell in output_cells:
      cell_type = cell_types_map[cell]
      if cell_type >= num_output_cell_types: continue
      cell_state = state[(cell * num_genes):((cell + 1) * num_genes)]
      assert len(cell_state) == num_genes
      output_states[i] += cell_state[-num_output_genes:]
  for state in output_states:
    assert len(state) == len(output_states[0])
  if all([ cell_type >= num_output_cell_types for cell_type in [ cell_types_map[cell] for cell in output_cells ] ]):
    print("\nNone of the output cells belong to the output cell type. Exiting...")
    sys.exit(1)

  x = output_states
  if function == 'parity':
    parity = [ 0 ] * (len(input_signal) + 1 - window_size)
    for i, j in enumerate(range(window_size, len(input_signal) + 1)):
      bitsum = 0
      for bit in input_signal[(j - window_size):j]:
        bitsum += bit
      if bitsum % 2 == 0:
        parity[i] = 0
      else:
        parity[i] = 1
    assert len(input_signal) == len(parity) + window_size - 1
    y = parity
  elif function == 'median':
    median = [ 0 ] * (len(input_signal) + 1 - window_size)
    for i, j in enumerate(range(window_size, len(input_signal) + 1)):
      bitsum = 0
      for bit in input_signal[(j - window_size):j]:
        bitsum += bit
      if bitsum > window_size // 2:
        median[i] = 1
      else:
        median[i] = 0
    assert len(input_signal) == len(median) + window_size - 1
    y = median

  if delay > 0:
    x = x[delay:]
    y = y[:-delay]

  if warmup_steps > 0:
    x = x[warmup_steps:]
    y = y[warmup_steps:]

  assert len(x) == len(y)
  for s in x:
    assert len(s) <= num_output_cells * num_output_genes

  return x, y, input_signal_info

def train_lasso(x, y, num_cells, num_output_cells, num_output_cell_types, num_output_genes, biocellion_output_file, output_dir, threads):
  """
    Train Lasso on the provided x and y lists and plots a histogram of features used from each cell and cell type.
  """
  from sklearn.linear_model import Lasso, LassoCV
  from sklearn.model_selection import train_test_split
  from sklearn.utils import parallel_backend

  x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.25)

  lasso = LassoCV(n_jobs=threads, selection='random', tol=0.05)
  lasso.fit(x_train, y_train)

  train_predicted = [ 1 if x > 0.5 else 0 for x in lasso.predict(x_train) ]
  test_predicted = [ 1 if x > 0.5 else 0 for x in lasso.predict(x_test) ]

  train_score = lasso.score(x_train, y_train)
  test_score = lasso.score(x_test, y_test)
  coeff_used = np.sum(lasso.coef_ != 0)
  coeffs = len(x[0])

  train_accuracy = sum([ 1 if a == b else 0 for a, b in zip(train_predicted, y_train) ]) / len(train_predicted)
  test_accuracy = sum([ 1 if a == b else 0 for a, b in zip(test_predicted, y_test) ]) / len(test_predicted)

  cell_types_map = get_cell_types(biocellion_output_file)
  cell_feature_count = {}
  cell_type_feature_count = {}
  for cell in range(num_cells - 1, num_cells - num_output_cells - 1, -1):
    cell_type = cell_types_map[cell]
    if cell_type >= num_output_cell_types: continue
    coeff_start = (num_cells - 1 - cell) * num_output_genes
    coeff_end = (num_cells - cell) * num_output_genes
    feature_count = np.sum(lasso.coef_[coeff_start:coeff_end] != 0)
    if feature_count > 0:
      if cell in cell_feature_count:
        cell_feature_count[cell] += feature_count
      else:
        cell_feature_count[cell] = feature_count
      if cell_type in cell_type_feature_count:
        cell_type_feature_count[cell_type] += feature_count
      else:
        cell_type_feature_count[cell_type] = feature_count
  
  plt.figure()
  plt.bar(range(len(cell_feature_count)), sorted(cell_feature_count.values(), reverse=True))
  xticks = [math.ceil(x) for x in plt.xticks()[0]]
  xticks_new = []
  for xtick in xticks:
    if xtick < 0: continue
    if xtick >= len(cell_feature_count):
      xticks_new.append(len(cell_feature_count) - 1)
      break
    else:
      xticks_new.append(xtick)
  xticks = xticks_new
  plt.xticks(xticks)
  plt.yticks([math.ceil(y) for y in plt.yticks()[0]])
  plt.savefig(os.path.join(output_dir, "cell_feature_count.png"))

  plt.figure()
  plt.bar(range(len(cell_type_feature_count)), sorted(cell_type_feature_count.values(), reverse=True))
  xticks = [math.ceil(x) for x in plt.xticks()[0]]
  xticks_new = []
  for xtick in xticks:
    if xtick < 0: continue
    if xtick >= len(cell_type_feature_count):
      xticks_new.append(len(cell_type_feature_count) - 1)
      break
    else:
      xticks_new.append(xtick)
  xticks = xticks_new
  plt.xticks(xticks)
  plt.yticks([math.ceil(y) for y in plt.yticks()[0]])
  plt.savefig(os.path.join(output_dir, "cell_type_feature_count.png"))

  return train_accuracy, test_accuracy, coeff_used, coeffs