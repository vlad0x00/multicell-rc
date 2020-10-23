"""
  List of simulation params.
"""

import argparse

# [ 0, inf ] integer params
def zeroplus_int(s):
  i = int(s)
  if i < 0: raise argparse.ArgumentTypeError(s + " is invalid, please input an integer >=0.")
  return i

# [ 1, inf ] integer params
def abovezero_int(s):
  i = int(s)
  if i <= 0: raise argparse.ArgumentTypeError(s + " is invalid, please input an integer >0.")
  return i

# [ 0.0, inf  ] float params
def zeroplus_float(s):
  f = float(s)
  if f < 0: raise argparse.ArgumentTypeError(s + " is invalid, please input a float >=0.")
  return f

# ( 0.0, inf  ] float params
def abovezero_float(s):
  f = float(s)
  if f <= 0: raise argparse.ArgumentTypeError(s + " is invalid, please input a float >0.")
  return f

# [ 0.0, 1.0 ] float params
def fraction_type(s):
  f = float(s)
  if f < 0 or f > 1: raise argparse.ArgumentTypeError(s + " is invalid, please input a number >=0 and <=1.")
  return f

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-g', '--genes', type=abovezero_int, default=100, help="Number of (internal) genes per cell.")
parser.add_argument('-c', '--cells', type=abovezero_int, default=12*12*12, help="Number of cells in the simulation.")
parser.add_argument('-p', '--cell-types', type=abovezero_int, default=10, help="Number of cell types in the simulation.")
parser.add_argument('-D', '--tissue-depth', type=zeroplus_int, default=3, help="Depth of tissue, in square cell layers, from the input substance source to the output cell layer. If zero, cells are arranged in a cube.")
parser.add_argument('-G', '--output-gene-fraction', type=fraction_type, default=1.0, help="Fraction of (internal) genes used for output.")
parser.add_argument('-C', '--output-cell-fraction', type=fraction_type, default=1.0/3.0, help="Fraction of cells used for output.")
parser.add_argument('-P', '--output-cell-type-fraction', type=fraction_type, default=1.0, help="Fraction of cell types used for output.")
parser.add_argument('-o', '--output-cells-random', action='store_true', help="Use random cells from anywhere in the tissue for LASSO training instead of only from the end.")
parser.add_argument('-k', '--in-degree', type=abovezero_int, default=2, help="Average node in-degree of gene network(s).")
parser.add_argument('-l', '--input-fraction', type=fraction_type, default=0.5, help="Fraction of nodes connected to the input signal.")
parser.add_argument('-I', '--input-signal-depth', type=zeroplus_int, default=1, help="Average node in-degree of gene network(s).")
parser.add_argument('-a', '--alpha-cytokines', type=zeroplus_float, default=1.0/45.0, help="Molecular decay rate of cytokines.")
parser.add_argument('-b', '--beta-cytokines', type=zeroplus_float, default=5.0, help="Grid diffusion coefficient of cytokines.")
parser.add_argument('-y', '--cytokines', type=zeroplus_int, default=5, help="Number of cytokines in the simulation.")
parser.add_argument('-L', '--secretion-low', type=zeroplus_float, default=1.0, help="Basal cytokine secretion when the gene is off.")
parser.add_argument('-H', '--secretion-high', type=zeroplus_float, default=5.0, help="Cytokine secretion when the gene is on.")
parser.add_argument('-t', '--cytokine-threshold', type=zeroplus_float, default=11.5, help="Cytokine threshold to turn a gene on. Value is relative to basal secretion.")
parser.add_argument('-r', '--cell-radius', type=abovezero_float, default=10.0, help="Cell radius.")
parser.add_argument('-v', '--voxel-length', type=abovezero_float, default=20.0, help="Simulation space voxel edge length.")
parser.add_argument('-T', '--timesteps', type=abovezero_int, default=1000, help="Number of simulation timesteps.")
parser.add_argument('-W', '--warmup-timesteps', type=zeroplus_int, default=100, help="Number of initial simulation timesteps to exclude from training.")
parser.add_argument('-d', '--delay', type=zeroplus_int, default=2, help="Step delay between input signal and output layer prediction.")
parser.add_argument('-f', '--function', choices=[ 'median', 'parity' ], default="parity", help="Function to learn")
parser.add_argument('-w', '--window-size', type=abovezero_int, default=3, help="Window size of predicted functions.")
parser.add_argument('-K', '--kappa', action='store_true', help="Enable kappa in Biocellion. Kappa reduces available space for molecules by the space occupied by cells.")
parser.add_argument('-S', '--summary', action='store_true', help="Enable output of summary files. Includes min, avg, and max values of cytokines and input signal throughout the simulation universe.")
parser.add_argument('-s', '--source-dist', type=abovezero_int, default=1, help="Distance between the signal producing wall and tissue.")
parser.add_argument('-R', '--reuse', action='store_true', help="Use previously generated initial gene network states, functions, and input signal. Otherwise generate new.")
parser.add_argument('-X', '--auxiliary', action='store_true', help="Generate auxiliary files (for debugging and/or visualization) -- dot and png files of the gene network model and a dump of all simulation states in text form.")
parser.add_argument('-O', '--output', default="output", help="Path of simulation output directory")
parser.add_argument('-j', '--threads', type=abovezero_int, default=2, help="Number of threads to use for the run.")

# Parse the arguments
def parse_args(args=None):
  global parser
  return parser.parse_args(args)