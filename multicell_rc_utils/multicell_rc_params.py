"""
  List of simulation params.
"""

import argparse
from functools import partial
from enum import Enum


class Function(Enum):
    PARITY = "parity"
    MEDIAN = "median"
    ALL_THREE_BIT_RECURSIVE = "256-recursive"
    ALL_THREE_BIT = "256-nonrecursive"


# Parse the arguments
def parse_args(args=None):
    def int_range(strval, minval=0, maxval=2**64):
        val = int(strval)
        if val < minval or val > maxval:
            raise argparse.ArgumentTypeError(
                f"{val} is invalid, please input an integer >{minval} and <{maxval}."
            )
        return val

    def float_range(strval, minval=0.0, maxval=float("inf")):
        val = float(strval)
        if val < minval or val > maxval:
            raise argparse.ArgumentTypeError(
                f"{val} is invalid, please input a float >={minval} and <={maxval}."
            )
        return val

    fraction_type = partial(float_range, minval=0.0, maxval=1.0)
    nonnegative_int = partial(int_range, minval=0)
    positive_int = partial(int_range, minval=1)
    nonnegative_float = partial(float_range, minval=0.0)
    oneplus_float = partial(float_range, minval=1.0)

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-g",
        "--genes",
        type=positive_int,
        default=100,
        help="Number of (reservoir) genes per cell.",
    )
    parser.add_argument(
        "-c",
        "--cells",
        type=positive_int,
        default=12 * 12 * 12,
        help="Number of cells in the simulation.",
    )
    parser.add_argument(
        "-p",
        "--strains",
        type=positive_int,
        default=10,
        help="Number of strains in the simulation.",
    )
    parser.add_argument(
        "-D",
        "--tissue-depth",
        type=nonnegative_int,
        default=3,
        help="Depth of tissue, in square cell layers, from the input source to the output cell layer. If zero, cells are arranged in a cube.",
    )
    parser.add_argument(
        "-G",
        "--output-gene-fraction",
        type=fraction_type,
        default=1.0,
        help="Fraction of (reservoir) genes used for output.",
    )
    parser.add_argument(
        "-C",
        "--output-cell-fraction",
        type=fraction_type,
        default=1.0 / 3.0,
        help="Fraction of cells used for output.",
    )
    parser.add_argument(
        "-P",
        "--output-strain-fraction",
        type=fraction_type,
        default=1.0,
        help="Fraction of strains used for output.",
    )
    parser.add_argument(
        "-o",
        "--output-cells-random",
        action="store_true",
        help="Use random cells from anywhere in the tissue for LASSO training instead of only from the end.",
    )
    parser.add_argument(
        "-k",
        "--in-degree",
        type=oneplus_float,
        default=2.0,
        help="Average node in-degree of gene network(s).",
    )
    parser.add_argument(
        "-l",
        "--input-fraction",
        type=fraction_type,
        default=0.5,
        help="Fraction of nodes connected to the input signal.",
    )
    parser.add_argument(
        "-I",
        "--input-signal-depth",
        type=nonnegative_int,
        default=1,
        help="How many layers of cells are reached by the input signal.",
    )
    parser.add_argument(
        "-a",
        "--alpha-esm",
        type=nonnegative_float,
        default=1.0 / 45.0,
        help="Molecular decay rate of ESMs.",
    )
    parser.add_argument(
        "-b",
        "--beta-esm",
        type=nonnegative_float,
        default=5.0,
        help="Grid diffusion coefficient of ESMs.",
    )
    parser.add_argument(
        "-y",
        "--esms",
        type=nonnegative_int,
        default=5,
        help="Number of ESMs in the simulation.",
    )
    parser.add_argument(
        "-L",
        "--secretion-low",
        type=nonnegative_float,
        default=1.0,
        help="Basal ESM secretion when the gene is off.",
    )
    parser.add_argument(
        "-H",
        "--secretion-high",
        type=nonnegative_float,
        default=5.0,
        help="ESM secretion when the gene is on.",
    )
    parser.add_argument(
        "-t",
        "--esm-threshold",
        type=nonnegative_float,
        default=11.5,
        help="ESM threshold to turn a gene on. Value is relative to basal secretion.",
    )
    parser.add_argument(
        "-r", "--cell-radius", type=oneplus_float, default=10.0, help="Cell radius."
    )
    parser.add_argument(
        "-v",
        "--voxel-length",
        type=oneplus_float,
        default=20.0,
        help="Simulation space voxel edge length.",
    )
    parser.add_argument(
        "-T",
        "--timesteps",
        type=positive_int,
        default=1000,
        help="Number of simulation timesteps.",
    )
    parser.add_argument(
        "-W",
        "--warmup-timesteps",
        type=nonnegative_int,
        default=100,
        help="Number of initial simulation timesteps to exclude from training.",
    )
    parser.add_argument(
        "-d",
        "--delay",
        type=nonnegative_int,
        default=2,
        help="Step delay between input signal and output layer prediction.",
    )
    parser.add_argument(
        "-f",
        "--function",
        choices=[
            Function.MEDIAN.value,
            Function.PARITY.value,
            Function.ALL_THREE_BIT_RECURSIVE.value,
            Function.ALL_THREE_BIT.value,
        ],
        default=Function.PARITY.value,
        help="Function to learn.",
    )
    parser.add_argument(
        "-w",
        "--window-size",
        type=positive_int,
        default=3,
        help="Window size of predicted functions.",
    )
    parser.add_argument(
        "-K",
        "--kappa",
        action="store_true",
        help="Enable kappa in Biocellion. Kappa reduces available space for molecules by the space occupied by cells.",
    )
    parser.add_argument(
        "-S",
        "--summary",
        action="store_true",
        help="Enable output of summary files. Includes min, avg, and max values of ESMs and input signal throughout the simulation universe.",
    )
    parser.add_argument(
        "-s",
        "--source-dist",
        type=positive_int,
        default=1,
        help="Distance between the signal producing wall and tissue.",
    )
    parser.add_argument(
        "-R",
        "--reuse",
        action="store_true",
        help="Use previously generated initial gene network states, functions, and input signal. Otherwise generate new.",
    )
    parser.add_argument(
        "-X",
        "--auxiliary",
        action="store_true",
        help="Generate auxiliary files (for debugging and/or visualization) -- dot and png files of the gene network model and a dump of all simulation states in text form.",
    )
    parser.add_argument(
        "-e",
        "--lasso-weights",
        action="store_true",
        help="If provided, write regression weights to standard output.",
    )
    parser.add_argument(
        "-O", "--output", default="output", help="Path of simulation output directory"
    )
    parser.add_argument(
        "-j",
        "--threads",
        type=positive_int,
        default=4,
        help="Number of threads to use for the run.",
    )

    parsed_args = parser.parse_args(args)

    if (
        parsed_args.function == Function.ALL_THREE_BIT
        or parsed_args.function == Function.ALL_THREE_BIT_RECURSIVE
    ):
        if parsed_args.window_size != 3:
            raise ValueError("Window size must be 3 for 256 functions.")

    return parsed_args
