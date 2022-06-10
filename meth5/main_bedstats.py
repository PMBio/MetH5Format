import argparse
from pathlib import Path
from typing import List, Dict

import numpy as np
from meth5 import MetH5File
from meth5.util import argtype_M5File, argtype_genomic_range

__description__ = "List chromosomes and chunks in meth5 files"


def set_arguments(sc_args: argparse.ArgumentParser):
    sc_args.add_argument(
        "-i",
        "--input_m5_file",
        type=argtype_M5File,
        required=True,
        help="List of MetH5 files",
    )
    
    sc_args.add_argument(
        "-d",
        "--data",
        type=str,
        choices=["meth_rate", "cpgs_covered", "num_calls"],
        default="methylation",
        nargs="+",
        required=True,
        help="What data columns to plot - multiple choices possible",
    )
    
    sc_args.add_argument(
        "-b",
        "--bed",
        type=argparse.FileType("rt"),
        required=True,
        help="Bed file listing regions to compute statistics on",
    )
    
    sc_args.add_argument(
        "-l",
        "--llr_threshold",
        type=float,
        required=False,
        default=2.0,
        help="Log-likelihood ratio threshold to consider a call",
    )


def get_statistic(vals, data, llr_threshold):
    if data == "meth_rate":
        rates, _ = vals.get_llr_site_rate(llr_threshold=llr_threshold)
        return np.nanmean(rates)
    
    if data == "cpgs_covered":
        _, sites = vals.get_llr_site_rate(llr_threshold=llr_threshold)
        return len(sites)
    
    if data == "num_calls":
        calls, _ = vals.get_llr_site_aggregate(lambda llrs: np.sum(np.abs(llrs) > llr_threshold))
        return np.nansum(calls)


def main(input_m5_file: Path, bed: Dict, chunk_size: int, data: str, llr_threshold: float):
    with MetH5File(input_m5_file, "r", chunk_size=chunk_size) as f:
        print("chrom\tstart\tend\t" + "\t".join(data))
        for i, region in enumerate(bed, 1):
            region = region.split("\t")
            if len(region) < 3:
                raise ValueError(f"Invalid bed file content in line {i}")
            try:
                region = {"chrom": region[0], "start": int(region[1]), "end": int(region[2])}
            except:
                raise ValueError(f"Invalid bed file content in line {i} - can't parse start and end coordinate")
            if region["chrom"] not in f.get_chromosomes():
                raise ValueError(f"Chromosome {region['chrom']} not found in file.")

            vals = f[region["chrom"]].get_values_in_range(region["start"], region["end"])
            
            out_line = f"{region['chrom']}\t{region['start']}\t{region['end']}\t"
            out_data = "\t".join([str(get_statistic(vals, d, llr_threshold)) for d in data])
            out_line = out_line + out_data
            print(out_line, flush=True)
