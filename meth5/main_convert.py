import argparse
from pathlib import Path
from typing import List, Union, Optional

import numpy as np
import pysam
import pandas as pd
from meth5 import MetH5File

__description__ = "Convert one or more BAM files with MM and ML tags to a single MetH5 file"

def p_to_llr(p, prior=0.5, eps=1e-16):
    """
    Converts the posterior probability p(a|x) into a log-likelihood ratio
    log(p(x|a)/p(x|~a)) given a prior p(a)
    """
    p = np.clip(p, eps, 1 - eps)  # Avoid numerical errors when p is 0 or 1
    return -np.log(prior * (1 - p) / (p * (1 - prior)))


def set_arguments(sc_args: argparse.ArgumentParser):
    sc_args.add_argument("-i", "--input_bam_files", type=str, required=True, nargs="+", help="Input BAM file(s)")
    sc_args.add_argument("-o", "--output_m5_file", type=str, required=True, help="Output MetH5 file")
    sc_args.add_argument(
        "-c", "--chromosomes", type=str, required=False, nargs="+", help="List of chromosomes to be included"
    )


class BufferedMetH5Writer:
    def __init__(self, mf: MetH5File, buffer_size: int = 1000000):
        self.mf = mf
        self.chrom_buf = [None] * buffer_size
        self.read_names_buf = [None] * buffer_size
        self.start_buf = [None] * buffer_size
        self.llr_buf = [None] * buffer_size
        self.buffer_size = buffer_size
        self.i = 0
    
    def __enter__(self):
        return self
    
    def flush(self):
        df = pd.DataFrame(
            data={
                "chromosome": self.chrom_buf[: self.i],
                "read_name": self.read_names_buf[: self.i],
                "start": self.start_buf[: self.i],
                "end": self.start_buf[: self.i],
                "log_lik_ratio": self.llr_buf[: self.i],
            }
        )
        self.mf.add_to_h5_file(df, postpone_sorting_until_close=True)
        self.i = 0
    
    def __exit__(self, *args, **kwargs):
        self.flush()
    
    def write(self, chrom: str, read_name: str, start: int, llr: float):
        if self.i >= self.buffer_size:
            self.flush()
        self.chrom_buf[self.i] = chrom
        self.read_names_buf[self.i] = read_name
        self.start_buf[self.i] = start
        self.llr_buf[self.i] = llr
        self.i += 1


def convert_one_bam(
    writer: BufferedMetH5Writer, infile_path: Union[str, Path], chromosomes: Optional[List[str]] = None
):
    bam = pysam.AlignmentFile(infile_path)
    all_chromosomes = [bam.get_reference_name(i) for i in range(bam.nreferences)]
    if chromosomes is None:
        chromosomes = all_chromosomes
    for chrom in chromosomes:
        if chrom not in all_chromosomes:
            continue
        for read in bam.fetch(chrom):
            q_to_r = [r for q, r in read.aligned_pairs if q is not None]
            for mod_type in read.modified_bases_forward:
                for q, p in read.modified_bases_forward[mod_type]:
                    r = q_to_r[q]
                    if r is None:
                        continue
                    llr = p_to_llr(p / 256)
                    writer.write(chrom, read.query_name, r, llr)


def main(
    input_bam_files: List[Union[str, Path]],
    output_m5_file: Union[str, Path],
    chunk_size: int,
    chromosomes: Optional[List[str]] = None,
):
    with MetH5File(output_m5_file, "w", chunk_size=chunk_size) as f:
        with BufferedMetH5Writer(f) as bw:
            for infile_path in input_bam_files:
                try:
                    print("Reading ", infile_path)
                    convert_one_bam(bw, infile_path, chromosomes=chromosomes)
                except:
                    print("WARNING! Could not read", infile_path, " - Skipping")
        
        with MetH5File(output_m5_file, "a", chunk_size=chunk_size) as f:
            print("Creating index")
            f.create_chunk_index()
