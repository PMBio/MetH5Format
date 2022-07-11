import argparse
from pathlib import Path
from typing import List, Dict
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
        choices=["methylation", "coverage"],
        default="methylation",
        required=False,
        help="What to plot in the data column",
    )
    
    sc_args.add_argument(
        "-r",
        "--region",
        type=argtype_genomic_range,
        required=False,
        default=None,
        help="Genomic region formatted as chr:start-end",
    )
    
    sc_args.add_argument(
        "-n",
        "--track_name",
        type=str,
        required=False,
        default=None,
        help="Track name to be displayed in viewer",
    )


def print_bedgraph(f: MetH5File, region: Dict, data: str):
    if region["chrom"] not in f.get_chromosomes():
        raise ValueError(f"Chromosome {region['chrom']} not found in file.")
    
    chrom = f[region["chrom"]]
    if region["end"] == -1:
        region["end"] = chrom.get_genomic_range()[1]
    if data == "methylation":
        rates, ranges = chrom.get_values_in_range(region["start"], region["end"]).get_llr_site_rate()

    elif data == "coverage":
        rates, ranges = (chrom.get_values_in_range(region["start"], region["end"]).get_llr_site_aggregate(
            lambda x: sum(1 for xi in x if abs(xi) > 2)))

    for bs, r in zip(rates, ranges):
        print(f"{region['chrom']} {r[0]} {r[1]} {bs}")

def print_header(data:str, track_name: str):
    if track_name is None:
        track_name = data
    if data == "methylation":
        print(f"track type=bedGraph name={track_name} description=center_label visibility=display_mode color=252,127,44 altColor=25,4,248 graphType=heatmap viewLimits=0:1 midRange=0.50:0.50 midColor=255,255,255")
    elif data == "coverage":
        print(f"track type=bedGraph name={track_name} description=center_label visibility=display_mode graphType=bar ")
    
def main(input_m5_file: Path, region: Dict, chunk_size: int, data: str, track_name: str):
    with MetH5File(input_m5_file, "r", chunk_size=chunk_size) as f:
        print_header(data, track_name)
        if region is not None:
            # Print a single region
            print_bedgraph(f, region, data)
            return
        
        # No region was specified, print everything
        for chrom in f.get_chromosomes():
            region = {"chrom": chrom, "start":0, "end":-1}
            print_bedgraph(f, region, data)
    
