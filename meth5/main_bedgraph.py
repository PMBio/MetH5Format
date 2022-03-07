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
        required=True,
        help="What to plot in the data column",
    )
    
    sc_args.add_argument(
        "-r",
        "--region",
        type=argtype_genomic_range,
        required=True,
        help="Genomic region formatted as chr:start-end",
    )


def main(input_m5_file: Path, region: Dict, chunk_size: int, data: str):
    with MetH5File(input_m5_file, "r", chunk_size=chunk_size) as f:
        if region["chrom"] not in f.get_chromosomes():
            raise ValueError(f"Chromosome {region['chrom']} not found in file.")
        
        if data == "methylation":
            rates, ranges = f[region["chrom"]].get_values_in_range(region["start"], region["end"]).get_llr_site_rate()
            print(
                "track type=bedGraph name=methylation description=center_label visibility=display_mode color=252,127,44 altColor=25,4,248 graphType=heatmap viewLimits=0:1 midRange=0.50:0.50 midColor=255,255,255"
            )
        elif data == "coverage":
            rates, ranges = (
                f[region["chrom"]]
                .get_values_in_range(region["start"], region["end"])
                .get_llr_site_aggregate(lambda x: sum(1 for xi in x if abs(xi) > 2))
            )
            print(
                "track type=bedGraph name=coverage description=center_label visibility=display_mode color=252,127,"
                "44 altColor=25,4,248 graphType=heatmap viewLimits=0:1 midRange=0.50:0.50 midColor=255,255,255"
            )
        
        for bs, r in zip(rates, ranges):
            print(f"{region['chrom']} {r[0]} {r[1]} {bs}")
