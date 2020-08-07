import os
import random
import csv
import shutil
import subprocess
import argparse
import struct

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment
from xml.dom import minidom

def zeroplus_int(s):
  i = int(s)
  if i < 0: raise argparse.ArgumentTypeError(s + " is invalid, please input an integer >=0.")
  return i

def abovezero_int(s):
  i = int(s)
  if i <= 0: raise argparse.ArgumentTypeError(s + " is invalid, please input an integer >0.")
  return i

def zeroplus_float(s):
  f = float(s)
  if f < 0: raise argparse.ArgumentTypeError(s + " is invalid, please input a float >=0.")
  return f

def abovezero_float(s):
  f = float(s)
  if f <= 0: raise argparse.ArgumentTypeError(s + " is invalid, please input a float >0.")
  return f

def fraction_type(s):
  f = float(s)
  if f < 0 or f > 1: raise argparse.ArgumentTypeError(s + " is invalid, please input a number >=0 and <=1.")
  return f

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-g', '--genes', type=abovezero_int, default=20, help="Number of (internal) genes per cell.")
parser.add_argument('-c', '--cells', type=abovezero_int, default=216, help="Number of cells in the simulation.")
parser.add_argument('-p', '--cell-types', type=abovezero_int, default=12, help="Number of cell types in the simulation.")
parser.add_argument('-u', '--output-gene-fraction', type=fraction_type, default=0.5, help="Fraction of (internal) genes used for output.")
parser.add_argument('-d', '--output-cell-fraction', type=fraction_type, default=0.5, help="Fraction of cells used for output.")
parser.add_argument('-k', '--degree', type=abovezero_int, default=2, help="Average node in-degree of gene network(s).")
parser.add_argument('-l', '--input-fraction', type=fraction_type, default=1.0, help="Fraction of nodes connected to the input signal.")
parser.add_argument('-f', '--function', choices=[ "median", "parity" ], default="parity", help="Function to learn")
parser.add_argument('-a', '--alpha', type=zeroplus_float, default=0.55, help="Molecular decay rate.")
parser.add_argument('-b', '--beta', type=zeroplus_float, default=5.0, help="Grid diffusion coefficient.")
parser.add_argument('-y', '--cytokines', type=zeroplus_int, default=2, help="Number of cytokines in the simulation.")
parser.add_argument('-o', '--secretion-low', type=zeroplus_float, default=0.0, help="Cytokine secretion when the gene is off.")
parser.add_argument('-i', '--secretion-high', type=zeroplus_float, default=40.0, help="Cytokine secretion when the gene is on.")
parser.add_argument('-t', '--cytokine-threshold', type=zeroplus_float, default=1.25, help="Cytokine threshold to turn a gene on.")
parser.add_argument('-r', '--cell-radius', type=abovezero_float, default=1.00, help="Cell radius.")
parser.add_argument('-x', '--grid-spacing', type=abovezero_float, default=2.3, help="Simulation space voxel length.")
parser.add_argument('-s', '--steps', type=abovezero_int, default=300, help="Number of simulation steps.")
parser.add_argument('-m', '--memory', type=zeroplus_int, default=0, help="Step delay between input signal and output layer prediction.")
parser.add_argument('-w', '--window-size', type=abovezero_int, default=5, help="Window size of predicted functions.")
parser.add_argument('-e', '--reuse', action='store_true', help="Use previously generated initial gene network states, functions, and input signal. Otherwise generate new.")
parser.add_argument('-z', '--visualize', action='store_true', help="Generate a dot and png file of the gene network model and every state of the simulation for visualization and debugging.")
parser.add_argument('-q', '--output', default="output", help="Path of simulation output directory")
parser.add_argument('-j', '--threads', type=abovezero_int, default=2, help="Number of threads to use for the run.")

def parse_args(args=None):
  global parser
  return parser.parse_args(args)

STATES_FILE = 'states'

class Node:

  def __init__(self, name):
    self.name = name

  def __str__(self):
    return self.name

  def __hash__(self):
    return self.name.__hash__()

def make_network_dot(num_genes, varf, dot_file):
  import networkx as nx

  graph = nx.DiGraph()

  for cell_type, variable_matrix in enumerate(varf):
    for gene, _ in enumerate(variable_matrix):
      node = Node(str(cell_type) + "_" + str(gene))
      graph.add_node(node)

  for cell_type, variable_matrix in enumerate(varf):
    for gene, variables in enumerate(variable_matrix):
      node = Node(str(cell_type) + "_" + str(gene))
      for v in variables:
        var_node = Node(str(cell_type) + "_" + str(v))
        graph.add_edge(var_node, node)

  pydot_graph = nx.drawing.nx_pydot.to_pydot(graph)
  pydot_graph.set_strict(False)
  pydot_graph.set_name("gene_networks")
  pydot_graph.write(dot_file, prog='dot')

  if shutil.which('dot') is not None:
    filename, file_extension = os.path.splitext(dot_file)
    with open(filename + ".png", 'w') as f:
      if num_genes > 100:
        subprocess.run([ 'sfdp', '-x', '-Goverlap=scale', '-T', 'png', dot_file ], stdout=f)
      else:
        subprocess.run([ 'dot', '-T', 'png', dot_file ], stdout=f)

def generate_gene_functions(num_cell_types, num_genes, connectivity, input_connections, num_cytokines, nv_file, varf_file, tt_file, dot_file):
  assert input_connections < num_genes

  import numpy as np

  nv = []
  varf = []
  tt = []

  for _ in range(num_cell_types):
    nv.append(np.zeros(num_genes, dtype=np.int32))

    total_edges = num_genes * connectivity - input_connections
    edges = []
    for _ in range(total_edges):
      while True:
        i, j = random.sample(range(1, num_genes), 2)
        if 1 + num_cytokines <= i < 1 + 2 * num_cytokines: continue
        if j < 1 + num_cytokines: continue
        if not (i, j) in edges:
          nv[-1][j] += 1
          edges.append((i, j))
          break
    for _ in range(input_connections):
      while True:
        j = random.randrange(1 + num_cytokines, num_genes)
        if not (0, j) in edges:
          nv[-1][j] += 1
          edges.append((0, j))
          break

    assert connectivity - 0.05 < sum(nv[-1]) / len(nv[-1]) < connectivity + 0.05

    varf.append([])
    for gene, n in enumerate(nv[-1]):
      varf[-1].append([])
      for _ in range(n):
        idx = -1
        for i, edge in enumerate(edges):
          if edge[1] == gene:
            varf[-1][-1].append(edge[0])
            idx = i
            break
        edges.pop(idx)

    tt.append([])
    for gene, n in enumerate(nv[-1]):
      if n > 0:
        tt[-1].append(np.random.randint(2, size=(2 ** n)))
      else:
        tt[-1].append([])

  with open(nv_file, 'w') as f:
    for cell_type_nv in nv:
      f.write(' '.join([str(x) for x in cell_type_nv]))
      f.write('\n')

  with open(varf_file, 'w') as f:
    for cell_type_varf in varf:
      for gene_varf in cell_type_varf:
        f.write(' '.join([str(x) for x in gene_varf]))
        f.write('\n')

  with open(tt_file, 'w') as f:
    for cell_type_tt in tt:
      for gene_tt in cell_type_tt:
        f.write(' '.join([str(x) for x in gene_tt]))
        f.write('\n')

  if dot_file != None:
    make_network_dot(num_genes, varf, dot_file)

def generate_input_signal(signal_len, signal_file):
  import numpy as np
  arr = np.random.randint(2, size=signal_len)
  arr[0] = 0 # Initial substance level is 0, otherwise grid phi initialization is non-trivial
  with open(signal_file, 'w') as f:
    f.write(' '.join([str(x) for x in arr]))

def generate_gene_initial_states(num_genes, num_cells, num_cytokines, input_signal_file, state_file):
  import numpy as np
  with open(input_signal_file, 'r') as f:
    input_signal = [ int(x) for x in f.readline().split() ]
  with open(state_file, 'w') as f:
    for _ in range(num_cells):
      state = np.random.randint(2, size=num_genes)
      state[0] = input_signal[0] # 0
      for i in range(num_cytokines):
        state[i + 1] = 0 # Initial substance levels are 0, otherwise grid phi initialization is non-trivial
      f.write(' '.join([str(x) for x in state]))
      f.write('\n')

def get_gene_values(output_file, num_genes, num_cells):
  with open(output_file, "rb") as f:
    lines = f.readlines()

  raw_start = 0
  raw_end = 0
  cell_id_offsets = []
  genebits_offsets = []
  genebits_updated = False
  for i, line in enumerate(lines):
    if b'<AppendedData encoding=\"raw\">' in line: raw_start = i + 1
    if b'</AppendedData>' in line: raw_end = i
    if b'offset' in line:
      tmp = line.decode('utf-8').split()
      for token in tmp:
        if 'offset' in token:
          tmp = token.split('=')
      assert tmp[0] == 'offset'
      tmp = tmp[1][1:]
      for i in range(len(tmp)):
        if not tmp[i].isdigit():
          offset = int(tmp[:i])
          break
      if b'Name="color"' in line:
        cell_id_offsets.append(offset)
      if b'Name="genebits_' in line:
        genebits_offsets.append(offset)
        genebits_updated = True
      elif genebits_updated:
        genebits_offsets.append(offset)
        genebits_updated = False

  raw_data = b''
  for line in lines[raw_start:raw_end]:
    raw_data += line
  raw_data = raw_data[1:] # Remove the underscore
  raw_data = raw_data[:-1] # Remove the newline

  if genebits_updated:
    genebits_offsets.append(len(raw_data))
    genebits_updated = False

  cell_id_offsets.append(genebits_offsets[0])

  cell_ids = []
  for i in range(num_cells - 1, -1, -1):
    id_start = cell_id_offsets[1] - i * 8 - 8
    id_end = cell_id_offsets[1] - i * 8
    cell_ids.append(int(struct.unpack('d', raw_data[id_start:id_end])[0]))

  gene_values = []
  for _ in range(num_cells):
    gene_values.append([])

  gene_num = 0
  for i in range(1, len(genebits_offsets)):
    data_end = genebits_offsets[i]
    data_start = data_end - num_cells * 8

    genebits = raw_data[data_start:data_end]
    for cell in range(num_cells):
      for bitPos in range(64):
        if gene_num + bitPos >= num_genes: break
        val = (genebits[(cell * 8):((cell + 1) * 8)][bitPos // 8] & (1 << (bitPos % 8))) >> (bitPos % 8)
        gene_values[cell].append(val)
    if gene_num + bitPos >= num_genes: break
    gene_num += 64

  reordered_gene_values = []
  for cell in range(num_cells):
    reordered_gene_values.append(gene_values[cell_ids.index(cell)])

  return reordered_gene_values

def get_states(num_genes, num_output_genes, num_cells, window_size, timesteps, output_dir):
  states = []
  for step in range(timesteps + 1):
    states.append([])
    gene_values = get_gene_values(os.path.join(output_dir, 'agent_0_0_0_' + str(step) + '.vtp'), num_genes, num_cells)
    for values in gene_values:
      states[-1] += values

  with open(os.path.join(output_dir, STATES_FILE), 'w') as f:
    for state in states:
      for i, bit in enumerate(state):
        if i > 0:
          f.write(' ')
        f.write(str(bit))
      f.write('\n')

  return states

def train_lasso(input_signal_file, biocellion_output_file, output_dir, num_genes, num_cells, num_output_genes, num_output_cells, window_size, delay, timesteps, function, visualize, threads):
  import numpy as np
  from sklearn.linear_model import Lasso, LassoCV
  from sklearn.model_selection import train_test_split
  from sklearn.utils import parallel_backend

  with open(input_signal_file) as f:
    input_signal = [ int(x) for x in f.readline().split() ]

  states = get_states(num_genes, num_output_genes, num_cells, window_size, timesteps, output_dir)

  cell_input_matches = []
  for _ in range(num_cells):
    cell_input_matches.append([])
  for signal, genes in zip(input_signal, states):
    for cell in range(num_cells):
      cell_input_matches[cell].append(signal == genes[cell * num_genes])
  cells_correct_input = 0
  for cell_matches in cell_input_matches:
    if all(cell_matches):
      cells_correct_input += 1
  assert cells_correct_input > 0

  states = states[window_size:]
  for state in states:
    assert len(state) == num_cells * num_genes

  def make_simulation_dots(states, output_dir):
    import networkx as nx

    class Node:

      def __init__(self, name):
        self.name = name

      def __str__(self):
        return self.name

      def __hash__(self):
        return self.name.__hash__()

    num_genes = len(states[0])
    dot_available = shutil.which('dot') is not None
    for idx, state in enumerate(states):
      graph = nx.DiGraph()
      for i, bit in enumerate(state):
        graph.add_node(Node(str(i)), label=bit)

      filename = os.path.join(output_dir, "state" + str(idx).zfill(len(str(abs(timesteps)))))

      pydot_graph = nx.drawing.nx_pydot.to_pydot(graph)
      pydot_graph.set_strict(False)
      pydot_graph.set_name("state" + str(idx))
      pydot_graph.write(filename + ".dot", prog='dot')

      if dot_available:
        with open(filename + ".png", 'w') as f:
          if num_genes > 100:
            subprocess.run([ 'sfdp', '-x', '-Goverlap=scale', '-T', 'png', filename + ".dot" ], stdout=f)
          else:
            subprocess.run([ 'dot', '-T', 'png', filename + ".dot" ], stdout=f)

  if visualize:
    make_simulation_dots(states, output_dir)
  else:
    for filename in os.listdir(output_dir):
      if filename.startswith("state") and filename.endswith('.dot'):
        os.remove(os.path.join(output_dir, filename))

  parity = []
  for i in range(window_size, len(input_signal)):
    bitsum = 0
    for bit in input_signal[(i - window_size):i]:
      bitsum += bit
    if bitsum % 2 == 0:
      parity.append(0)
    else:
      parity.append(1)
  assert len(input_signal) == len(parity) + window_size

  median = []
  for i in range(window_size, len(input_signal)):
    bitsum = 0
    for bit in input_signal[(i - window_size):i]:
      bitsum += bit
    if bitsum > window_size // 2:
      median.append(1)
    else:
      median.append(0)
  assert len(input_signal) == len(median) + window_size

  functions = { 'parity' : parity, 'median' : median }

  output_states = []
  for state in states:
    output_states.append([])
    for cell in range(num_output_cells):
      if cell == 0:
        cell_state = state[-num_genes:]
      else:
        cell_state = state[-((cell + 1) * num_genes):-(cell * num_genes)]
      assert len(cell_state) == num_genes
      output_states[-1] += cell_state[-num_output_genes:]
  for state in output_states:
    assert len(state) == len(output_states[0])

  x = output_states
  y = functions[function]

  if delay > 0:
    x = x[delay:]
    y = y[:-delay]

  assert len(x) == len(y)
  for s in x:
    assert len(s) == num_output_cells * num_output_genes

  x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.25)

  lasso = LassoCV(n_jobs=threads)
  #lasso = Lasso(alpha=LASSO_ALPHA)
  lasso.fit(x_train, y_train)

  train_predicted = [ 1 if x > 0.5 else 0 for x in lasso.predict(x_train) ]
  test_predicted = [ 1 if x > 0.5 else 0 for x in lasso.predict(x_test) ]

  train_score = lasso.score(x_train, y_train)
  test_score = lasso.score(x_test, y_test)
  coeff_used = np.sum(lasso.coef_!=0)

  train_accuracy = sum([ 1 if a == b else 0 for a, b in zip(train_predicted, y_train) ]) / len(train_predicted)
  test_accuracy = sum([ 1 if a == b else 0 for a, b in zip(test_predicted, y_test) ]) / len(test_predicted)

  return train_accuracy, test_accuracy, coeff_used, cells_correct_input

def prettify(elem):
  """Return a pretty-printed XML string for the Element.
  """
  rough_string = ET.tostring(elem, 'utf-8')
  reparsed = minidom.parseString(rough_string)
  return reparsed.toprettyxml(indent="  ")

def make_params_xml(xml_path, output_dir, simulation_steps, additional_params, threads):
  xml_file = open(xml_path, 'w')

  # Biocellion required parameters
  bcell_num_baseline = simulation_steps
  bcell_nx = '10'
  bcell_ny = '10'
  bcell_nz = '10'
  bcell_partition_size = 10
  bcell_path = output_dir
  bcell_interval = 1
  bcell_start_x = 0
  bcell_start_y = 0
  bcell_start_z = 0
  bcell_size_x = 10
  bcell_size_y = 10
  bcell_size_z = 10

  # Biocellion optional parameteres
  bcell_input_param = additional_params
  bcell_verbosity = 0 # [0-5]

  bcell_num_threads = threads
  bcell_num_node_groups = 1
  bcell_num_nodes_per_group = 1
  bcell_num_sockets_per_node = 1
  bcell_max_load_imbalance = 1.2

  bcell_super_x = bcell_size_x
  bcell_super_y = bcell_size_y
  bcell_super_z = bcell_size_z

  bcell_summary = 1
  bcell_load_balance = 0
  bcell_regridding = 1000
  bcell_checkpoint = 0

  bcell_refine_ratio = 2
  bcell_fill_ratio = 0.5

  #  Write xml for Biocellion
  top = Element('biocellion')

  comment = Comment('Input parameters for Biocellion, generated by Python parser.')
  top.append(comment)

  # Printing required parameters
  xml_required = SubElement(top, 'required')
  xml_required.append(Comment('Required Parameters for biocellion.'))

  xml_time_step = SubElement(xml_required, 'time_step')
  xml_time_step.set('num_baseline', str(bcell_num_baseline))

  xml_domain = SubElement(xml_required, 'domain')
  xml_domain.set('x', str(bcell_nx))
  xml_domain.set('y', str(bcell_ny))
  xml_domain.set('z', str(bcell_nz))

  xml_init_data = SubElement(xml_required, 'init_data')
  xml_init_data.set('partition_size', str(bcell_partition_size))
  xml_init_data.set('src', 'code')

  xml_output = SubElement(xml_required, 'output ')
  xml_output.set('path', bcell_path)
  xml_output.set('interval', str(bcell_interval))
  xml_output.set('particle', 'pvtp')
  xml_output.set('grid', 'vthb')
  xml_output.set('start_x', str(bcell_start_x))
  xml_output.set('start_y', str(bcell_start_y))
  xml_output.set('start_z', str(bcell_start_z))
  xml_output.set('size_x', str(bcell_size_x))
  xml_output.set('size_y', str(bcell_size_y))
  xml_output.set('size_z', str(bcell_size_z))

  # Printing Required parameterss
  xml_optional = SubElement(top, 'optional')
  xml_optional.append(Comment('Optional parameters for biocellion.'))

  xml_bioparm_model = SubElement(xml_optional, 'model')
  xml_bioparm_model.set('param', bcell_input_param)

  xml_bioparm_stdout = SubElement(xml_optional, 'stdout')
  xml_bioparm_stdout.set('verbosity', str(bcell_verbosity))

  xml_bioparm_system = SubElement(xml_optional, 'system')
  xml_bioparm_system.set('num_node_groups', str(bcell_num_node_groups))
  xml_bioparm_system.set('num_nodes_per_group', str(bcell_num_nodes_per_group))
  xml_bioparm_system.set('num_sockets_per_node', str(bcell_num_sockets_per_node))
  xml_bioparm_system.set('max_load_imbalance', str(bcell_max_load_imbalance))#
  xml_bioparm_system.set('num_threads', str(bcell_num_threads))

  xml_bioparm_super_partition = SubElement(xml_optional, 'super_partition')
  xml_bioparm_super_partition.set('x', str(bcell_super_x))
  xml_bioparm_super_partition.set('y', str(bcell_super_y))
  xml_bioparm_super_partition.set('z', str(bcell_super_z))

  xml_interval = SubElement(xml_optional, 'interval')
  xml_interval.set('summary', str(bcell_summary))
  xml_interval.set('load_balance', str(bcell_load_balance))
  xml_interval.set('regridding', str(bcell_regridding))
  xml_interval.set('checkpoint', str(bcell_checkpoint))

  xml_amr = SubElement(xml_optional, 'amr')
  xml_amr.set('refine_ratio', str(bcell_refine_ratio))
  xml_amr.set('fill_ratio', str(bcell_fill_ratio))

  xml_file.write(prettify(top))
