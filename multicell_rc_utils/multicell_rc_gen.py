import os
import sys
import random
import shutil
import subprocess
import numpy as np
import networkx as nx
import math

import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment
from xml.dom import minidom

class Node:

  def __init__(self, name):
    self.name = name

  def __str__(self):
    return self.name

  def __hash__(self):
    return self.name.__hash__()

def make_network_dot(num_genes, varf, dot_file):
  """
    Generate the gene network .dot file for visualization.
  """
  graph = nx.DiGraph()

  # Generate the graph nodes
  for strain, variable_matrix in enumerate(varf):
    for gene, _ in enumerate(variable_matrix):
      node = Node(str(strain) + "_" + str(gene))
      graph.add_node(node)

  # Generate the graph edges
  for strain, variable_matrix in enumerate(varf):
    for gene, variables in enumerate(variable_matrix):
      node = Node(str(strain) + "_" + str(gene))
      for v in variables:
        var_node = Node(str(strain) + "_" + str(v))
        graph.add_edge(var_node, node)

  # Save the .dot file
  pydot_graph = nx.drawing.nx_pydot.to_pydot(graph)
  pydot_graph.set_strict(False)
  pydot_graph.set_name("gene_networks")
  pydot_graph.write(dot_file, prog='dot')

  # If the `dot` program exists, generate the graph image
  if shutil.which('dot') is not None:
    filename, file_extension = os.path.splitext(dot_file)
    with open(filename + ".png", 'w') as f:
      if num_genes > 100: # Large graphs are generated faster with sfdp
        subprocess.run([ 'sfdp', '-x', '-Goverlap=scale', '-T', 'png', dot_file ], stdout=f)
      else:
        subprocess.run([ 'dot', '-T', 'png', dot_file ], stdout=f)

def generate_ESM_normalization_files(num_genes, nv_file, varf_file, tt_file, signal_file, signal_len, state_file):
  nv = [ np.zeros(num_genes, dtype=np.int64) ]
  varf = [ [ [ 0 ] for _ in range(num_genes) ] ]
  tt = [ [ [ 0 ] for _ in range(num_genes) ] ]

  # Save nv
  with open(nv_file, 'w') as f:
    for strain_nv in nv:
      f.write(' '.join([str(x) for x in strain_nv]))
      f.write('\n')

  # Save varf
  with open(varf_file, 'w') as f:
    for strain_varf in varf:
      for gene_varf in strain_varf:
        f.write(' '.join([str(x) for x in gene_varf]))
        f.write('\n')

  # Save nv
  with open(tt_file, 'w') as f:
    for strain_tt in tt:
      for gene_tt in strain_tt:
        f.write(' '.join([str(x) for x in gene_tt]))
        f.write('\n')

  # Save input signal.
  arr = np.ones(signal_len, dtype=np.int64)
  with open(signal_file, 'w') as f:
    f.write(' '.join([str(x) for x in arr]))

  # Save initial states
  with open(state_file, 'w') as f:
    state = np.zeros(num_genes, dtype=np.int64)
    f.write(' '.join([str(x) for x in state]))
    f.write('\n')

def generate_gene_functions(num_strains, num_genes, connectivity, input_connections, num_ESMs, nv_file, varf_file, tt_file, dot_file):
  """
    Generate the truth tables for genes for each strain
  """
  assert input_connections < num_genes

  """
    Each of the three lists is first indexed by strain and then by gene.

    nv : The number of input gene variables for each gene function
    varf : The list of input genes variables for each gene function
    tt : The truth tables for each gene function
  """
  nv = []
  varf = []
  tt = []

  for strain in range(num_strains):
    # Initialize nv-s for current strain
    nv.append(np.zeros(num_genes, dtype=np.int32))

    # Ensure the number of edges corresponds to the average input degree of nodes
    total_edges = math.ceil((num_genes - 1 - num_ESMs) * connectivity)

    # Number of possible edges excluding the input gene
    possible_edges_no_input = (num_genes - 1 - num_ESMs) ** 2
    edges_no_input = [ (0, 0) ] * possible_edges_no_input
    idx = 0
    for i in range(1, num_genes):
      for j in range(1, num_genes):
        if j < (1 + num_ESMs): continue
        if (1 + num_ESMs) <= i < (1 + 2 * num_ESMs): continue
        edges_no_input[idx] = (i, j)
        idx += 1
    assert idx == possible_edges_no_input

    # Number of possible edges from the input gene
    possible_edges_input = num_genes - 1 - num_ESMs
    edges_input = [ (0, 0) ] * possible_edges_input
    idx = 0
    for j in range(1 + num_ESMs, num_genes):
      edges_input[idx] = (0, j)
      idx += 1
    assert idx == possible_edges_input

    # Shuffle the possible edges for a randomized network
    random.shuffle(edges_no_input)
    random.shuffle(edges_input)

    # And take a number of them
    edges = edges_no_input[:(total_edges - input_connections)]
    edges += edges_input[:input_connections]

    # Update nv
    for edge in edges:
      nv[-1][edge[1]] += 1

    # Ensure the average in-degree condition is satisfied
    assert (connectivity - 0.01) < (sum(nv[-1]) / (len(nv[-1]) - num_ESMs - 1)) < (connectivity + 0.01)

    # Update gene variables
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

    # Update truth tables
    tt.append([])
    for gene, n in enumerate(nv[-1]):
      if n > 0:
        tt[-1].append(np.random.randint(2, size=(2 ** n)))
      else:
        tt[-1].append([])

  # Save nv
  with open(nv_file, 'w') as f:
    for strain_nv in nv:
      f.write(' '.join([str(x) for x in strain_nv]))
      f.write('\n')

  # Save varf
  with open(varf_file, 'w') as f:
    for strain_varf in varf:
      for gene_varf in strain_varf:
        f.write(' '.join([str(x) for x in gene_varf]))
        f.write('\n')

  # Save nv
  with open(tt_file, 'w') as f:
    for strain_tt in tt:
      for gene_tt in strain_tt:
        f.write(' '.join([str(x) for x in gene_tt]))
        f.write('\n')

  # If dot_file is not None, then a .dot file should be saved
  if dot_file != None:
    make_network_dot(num_genes, varf, dot_file)

def generate_input_signal(signal_len, signal_file):
  """
    Generate a random binary signal.
  """
  arr = np.random.randint(2, size=signal_len)
  arr[0] = 0 # Initial substance level is 0, otherwise grid phi initialization is non-trivial
  with open(signal_file, 'w') as f:
    f.write(' '.join([str(x) for x in arr]))

def generate_gene_initial_states(num_genes, num_cells, num_ESMs, input_signal_file, state_file):
  """
    Generate random initial states for genes for all cells
  """
  # Input signal gene should be initialized with the same value as the signal
  with open(input_signal_file, 'r') as f:
    input_signal = [ int(x) for x in f.readline().split() ]
  with open(state_file, 'w') as f:
    for _ in range(num_cells):
      state = np.random.randint(2, size=num_genes)
      state[0] = input_signal[0] # 0
      for i in range(num_ESMs):
        state[i + 1] = 0 # Initial substance levels are 0, otherwise grid phi initialization is non-trivial
      f.write(' '.join([str(x) for x in state]))
      f.write('\n')

def prettify(elem):
  """
    Return a pretty-printed XML string for the Element.
  """
  rough_string = ET.tostring(elem, 'utf-8')
  reparsed = minidom.parseString(rough_string)
  return reparsed.toprettyxml(indent="  ")

def make_params_xml(xml_path, output_dir, simulation_steps, additional_params, threads, universe_length, universe_width):
  xml_file = open(xml_path, 'w')

  # Biocellion required parameters
  bcell_num_baseline = simulation_steps
  bcell_nx = str(max(4, universe_width))
  bcell_ny = str(max(4, universe_width))
  bcell_nz = str(max(4, universe_length))
  bcell_partition_size = max(4, universe_width, universe_length)
  bcell_path = output_dir
  bcell_interval = 1
  bcell_start_x = 0
  bcell_start_y = 0
  bcell_start_z = 0
  bcell_size_x = max(4, universe_width)
  bcell_size_y = max(4, universe_width)
  bcell_size_z = max(4, universe_length)

  # Biocellion optional parameteres
  bcell_input_param = additional_params
  bcell_verbosity = 0 # [0-5]

  bcell_num_threads = threads
  bcell_num_node_groups = 1
  bcell_num_nodes_per_group = 1
  bcell_num_sockets_per_node = 1
  bcell_max_load_imbalance = 1.2

  supersize = max(bcell_size_x, bcell_size_y, bcell_size_z)
  bcell_super_x = supersize
  bcell_super_y = supersize
  bcell_super_z = supersize

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
