import argparse
from pathlib import Path
from typing import List
from meth5 import MetH5File
from meth5.util import argtype_M5File

__description__ = "List chromosomes and chunks in meth5 files"


def set_arguments(sc_args: argparse.ArgumentParser):
    sc_args.add_argument(
        "-i",
        "--input_m5_files",
        type=argtype_M5File,
        required=True,
        nargs="+",
        help="List of MetH5 files",
    )


def pad(string: str, width: int, pad_char: str = " ", align: str = "left") -> str:
    if len(string) >= width:
        return string
    else:
        if align == "left":
            return string + pad_char * (width - len(string))
        elif align == "center":
            npad = width - len(string)
            if npad % 2 == 0:
                npad = npad + 1
            return (pad_char * (npad // 2)) + string + (pad_char * (npad // 2 + 1))
        elif align == "right":
            return pad_char * (width - len(string)) + string


def main(input_m5_files: List[Path], chunk_size: int):
    for m5file in input_m5_files:
        rows = [[" Chromosome ", " Number of chunks "]]
        
        with MetH5File(m5file, "r", chunk_size=chunk_size) as f:
            for chrom in f.get_chromosomes():
                rows.append([f" {chrom}", f" {f[chrom].get_number_of_chunks()}"])
        
        colwidths = [max((len(r[i]) for r in rows)) for i in range(2)]
        rows = [[pad(row[i], colwidths[i]) for i in range(len(row))] for row in rows]
        print(pad(f" {m5file.name} ", sum(colwidths) + len(colwidths), "-", align="center"))
        for row in rows:
            print(f"|{'|'.join(row)}|")
    print(pad("", sum(colwidths) + len(colwidths) + 1, "-", align="center"))
