import argparse

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
parser.add_argument('-g', '--genes', type=abovezero_int, default=100, help="Number of (internal) genes per cell.")
parser.add_argument('-c', '--cells', type=abovezero_int, default=3600, help="Number of cells in the simulation.")
parser.add_argument('-p', '--cell-types', type=abovezero_int, default=1000, help="Number of cell types in the simulation.")
parser.add_argument('-E', '--tissue-depth', type=abovezero_int, default=3, help="Depth of tissue, in square cell layers, from the input substance source to the output cell layer.")
parser.add_argument('-G', '--output-gene-fraction', type=fraction_type, default=0.5, help="Fraction of (internal) genes used for output.")
parser.add_argument('-C', '--output-cell-fraction', type=fraction_type, default=0.33333333, help="Fraction of cells used for output.")
parser.add_argument('-P', '--output-cell-type-fraction', type=fraction_type, default=0.5, help="Fraction of cell types used for output.")
parser.add_argument('-k', '--degree', type=abovezero_int, default=2, help="Average node in-degree of gene network(s).")
parser.add_argument('-l', '--input-fraction', type=fraction_type, default=1.0, help="Fraction of nodes connected to the input signal.")
parser.add_argument('-D', '--dirichlet-boundary', type=abovezero_float, default=5.0, help="Value of dirichlet boundary when input signal is on. The value is 0 when the signal is off.")
parser.add_argument('-T', '--input-threshold', type=abovezero_float, default=1.0, help="Threshold for the molecular concentration of input signal for cells to consider it on.")
parser.add_argument('-A', '--alpha-input', type=zeroplus_float, default=0.55, help="Molecular decay rate of input signal.")
parser.add_argument('-B', '--beta-input', type=zeroplus_float, default=5.0, help="Grid diffusion coefficient of input signal.")
parser.add_argument('-f', '--function', choices=[ "median", "parity" ], default="parity", help="Function to learn")
parser.add_argument('-a', '--alpha-cytokines', type=zeroplus_float, default=0.55, help="Molecular decay rate of cytokines.")
parser.add_argument('-b', '--beta-cytokines', type=zeroplus_float, default=5.0, help="Grid diffusion coefficient of cytokines.")
parser.add_argument('-y', '--cytokines', type=zeroplus_int, default=8, help="Number of cytokines in the simulation.")
parser.add_argument('-L', '--secretion-low', type=zeroplus_float, default=0.0, help="Cytokine secretion when the gene is off.")
parser.add_argument('-H', '--secretion-high', type=zeroplus_float, default=55.0, help="Cytokine secretion when the gene is on.")
parser.add_argument('-t', '--cytokine-threshold', type=zeroplus_float, default=1.5, help="Cytokine threshold to turn a gene on.")
parser.add_argument('-r', '--cell-radius', type=abovezero_float, default=1.00, help="Cell radius.")
parser.add_argument('-x', '--grid-spacing', type=abovezero_float, default=2.3, help="Simulation space voxel edge length.")
parser.add_argument('-s', '--steps', type=abovezero_int, default=1000, help="Number of simulation steps.")
parser.add_argument('-W', '--warmup-steps', type=zeroplus_int, default=100, help="Number of initial simulation steps to exclude from training.")
parser.add_argument('-m', '--memory', type=zeroplus_int, default=2, help="Step delay between input signal and output layer prediction.")
parser.add_argument('-w', '--window-size', type=abovezero_int, default=5, help="Window size of predicted functions.")
parser.add_argument('-R', '--reuse', action='store_true', help="Use previously generated initial gene network states, functions, and input signal. Otherwise generate new.")
parser.add_argument('-X', '--auxiliary', action='store_true', help="Generate auxiliary files (for debugging and/or visualization) -- dot and png files of the gene network model and a dump of all simulation states in text form.")
parser.add_argument('-O', '--output', default="output", help="Path of simulation output directory")
parser.add_argument('-j', '--threads', type=abovezero_int, default=2, help="Number of threads to use for the run.")

def parse_args(args=None):
  global parser
  return parser.parse_args(args)