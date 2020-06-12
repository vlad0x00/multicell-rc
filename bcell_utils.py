import numpy as np
import random
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element, SubElement, Comment
from xml.dom import minidom
import csv

def generate_gene_functions(nv_file, varf_file, tt_file, num_genes):
  nv = np.random.randint(num_genes, size=num_genes)

  varf = np.ones((num_genes, num_genes), dtype=np.int32) * -1
  for i, n in enumerate(nv):
    varf[i, :n] = random.sample(range(0, num_genes), n)

  tt = np.ones((num_genes, 2 ** num_genes), dtype=np.int32) * -1
  for i, n in enumerate(nv):
    if n == 0: continue
    tt[i, :(2 ** n)] = np.random.randint(2, size=(2 ** n))

  with open(nv_file, 'w') as f:
    f.write(' '.join([str(x) for x in nv]))

  with open(varf_file, 'w') as f:
    for row in varf:
      f.write(' '.join([str(x).rjust(2) for x in row]) + '\n')

  with open(tt_file, 'w') as f:
    for row in tt:
      f.write(' '.join([str(x).rjust(2) for x in row]) + '\n')

def generate_bn_initial_state(state_file, num_genes):
  state = np.random.randint(2, size=num_genes)
  with open(state_file, 'w') as f:
    f.write(' '.join([str(x) for x in state]))

def generate_input_array(array_file, array_len):
  arr = np.random.randint(2, size=array_len)
  with open(array_file, 'w') as f:
    f.write(' '.join([str(x) for x in arr]))

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
  bcell_nx = '16'
  bcell_ny = '16'
  bcell_nz = '16'
  bcell_partition_size = 16
  bcell_path = output_dir
  bcell_interval = 1
  bcell_start_x = 0
  bcell_start_y = 0
  bcell_start_z = 0
  bcell_size_x = 16
  bcell_size_y = 16
  bcell_size_z = 16

  # Biocellion optional parameteres
  bcell_input_param = additional_params
  bcell_verbosity = 0 # [0-5]

  bcell_num_threads = 6
  bcell_num_node_groups = 1
  bcell_num_nodes_per_group = 1
  bcell_num_sockets_per_node = 1
  bcell_max_load_imbalance = 1.2

  bcell_super_x = 64
  bcell_super_y = 64
  bcell_super_z = 64

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