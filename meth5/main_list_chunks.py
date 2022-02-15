import argparse
from pathlib import Path
from typing import List
from meth5 import MetH5File
from meth5.util import argtype_M5File

__description__ = "List chromosomes and chunks in meth5 files"

def set_arguments(sc_args: argparse.ArgumentParser):
    sc_args.add_argument(
        "--input_m5_files",
        type=argtype_M5File,
        required=True,
        nargs="+",
        help="List of MetH5 files",
    )

def main(m5files: List[Path], chunk_size: int):
    for m5file in m5files:
        print(f"{m5file.basename()}: ")
        with MetH5File(m5file, "r", chunk_size=chunk_size) as f:
            for chrom in f.get_chromosomes():
                print(f"{chrom}: {f[chrom].get_number_of_chunks()}")
