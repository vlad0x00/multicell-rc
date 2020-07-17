import os
import numpy as np
import random
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment
from xml.dom import minidom
import csv
import networkx as nx
from sklearn.linear_model import Lasso, LassoCV
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np

class Node:

  def __init__(self, name):
    self.name = name

  def __str__(self):
    return self.name

  def __hash__(self):
    return self.name.__hash__()

def make_network_dot(varf, dot_file):
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
        graph.add_edge(node, var_node)

  pydot_graph = nx.drawing.nx_pydot.to_pydot(graph)
  pydot_graph.set_strict(False)
  pydot_graph.set_name("gene_networks")
  pydot_graph.write(dot_file, prog='dot')

def generate_gene_functions(num_cell_types, num_genes, connectivity, input_connections, num_cytokines, nv_file, varf_file, tt_file, dot_file):
  assert input_connections < num_genes

  nv = []
  varf = []
  tt = []

  for _ in range(num_cell_types):
    nv.append(np.zeros(num_genes, dtype=np.int32))

    total_edges = (num_genes - 1 - num_cytokines) * connectivity
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

    assert connectivity - 1 < sum(nv[-1]) / len(nv[-1]) < connectivity + 1

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
    make_network_dot(varf, dot_file)

def generate_input_signal(signal_len, signal_file):
  arr = np.random.randint(2, size=signal_len)
  with open(signal_file, 'w') as f:
    f.write(' '.join([str(x) for x in arr]))

def generate_gene_initial_states(num_genes, num_cells, input_signal_file, state_file):
  with open(input_signal_file, 'r') as f:
    input_signal = [ int(x) for x in f.readline().split() ]
  with open(state_file, 'w') as f:
    for _ in range(num_cells):
      state = np.random.randint(2, size=num_genes)
      state[0] = input_signal[0]
      f.write(' '.join([str(x) for x in state]))
      f.write('\n')

def get_gene_values(output_file, num_genes, num_cells):
  with open(output_file, "rb") as f:
    lines = f.readlines()

  raw_start = 0
  raw_end = 0
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

  if genebits_updated:
    genebits_offsets.append(len(raw_data))
    genebits_updated = False

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
        gene_values[cell].append((genebits[(cell * 8):((cell + 1) * 8)][bitPos // 8] & (1 << (bitPos % 8))) >> (bitPos % 8))
    if gene_num + bitPos >= num_genes: break
    gene_num += 64
  gene_values.reverse()

  return gene_values

def train_lasso(input_signal_file, biocellion_output_file, output_dir, num_genes, num_output_genes, num_cells, window_size, delay, timesteps, visualize):
  with open(input_signal_file) as f:
    input_signal = [ int(x) for x in f.readline().split() ]

  states = []
  for step in range(timesteps + 1):
    states.append([])
    gene_values = get_gene_values(os.path.join(output_dir, 'agent_0_0_0_' + str(step) + '.vtp'), num_genes, num_cells)
    for values in gene_values:
      states[-1] += values

  states = states[window_size:]
  for state in states:
    assert len(state) == num_cells * num_output_genes

  def make_simulation_dots(states, output_dir):
    class Node:

      def __init__(self, name):
        self.name = name

      def __str__(self):
        return self.name

      def __hash__(self):
        return self.name.__hash__()

    for idx, state in enumerate(states):
      graph = nx.DiGraph()
      for i, bit in enumerate(state):
        graph.add_node(Node(str(i)), label=bit)

      pydot_graph = nx.drawing.nx_pydot.to_pydot(graph)
      pydot_graph.set_strict(False)
      pydot_graph.set_name("state" + str(idx))
      pydot_graph.write(output_dir + "state" + str(index).zfill(3) + ".dot", prog='dot')
  if visualize:
    make_simulation_dots(states, output_dir)
  else:
    for filename in os.listdir(output_dir):
      if filename.startswith("state"):
        os.remove(output_dir + filename)

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

  x = states
  y = median

  if delay > 0:
    x = x[delay:]
    y = y[:-delay]

  assert len(x) == len(y)
  for s in x:
    assert len(s) == num_cells * num_genes

  x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.25)

  lasso = LassoCV(max_iter=10000)
  #lasso = Lasso(alpha=LASSO_ALPHA)
  lasso.fit(x_train, y_train)

  train_predicted = [ 1 if x > 0.5 else 0 for x in lasso.predict(x_train) ]
  test_predicted = [ 1 if x > 0.5 else 0 for x in lasso.predict(x_test) ]

  train_score = lasso.score(x_train, y_train)
  test_score = lasso.score(x_test, y_test)
  coeff_used = np.sum(lasso.coef_!=0)

  train_accuracy = sum([ 1 if a == b else 0 for a, b in zip(train_predicted, y_train) ]) / len(train_predicted)
  test_accuracy = sum([ 1 if a == b else 0 for a, b in zip(test_predicted, y_test) ]) / len(test_predicted)

  print("Training accuracy:", train_accuracy)
  print("Testing accuracy:", test_accuracy)
  print("Number of features used:", coeff_used)

def prettify(elem):
  """Return a pretty-printed XML string for the Element.
  """
  rough_string = ET.tostring(elem, 'utf-8')
  reparsed = minidom.parseString(rough_string)
  return reparsed.toprettyxml(indent="  ")

def make_params_xml(xml_path, output_dir, simulation_steps, additional_params):
  xml_file = open(xml_path, 'w')

  # Biocellion required parameters
  bcell_num_baseline = simulation_steps
  bcell_nx = '8'
  bcell_ny = '8'
  bcell_nz = '8'
  bcell_partition_size = 8
  bcell_path = output_dir
  bcell_interval = 1
  bcell_start_x = 0
  bcell_start_y = 0
  bcell_start_z = 0
  bcell_size_x = 8
  bcell_size_y = 8
  bcell_size_z = 8

  # Biocellion optional parameteres
  bcell_input_param = additional_params
  bcell_verbosity = 0 # [0-5]

  bcell_num_threads = 4
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