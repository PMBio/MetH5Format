from pathlib import Path
import argparse
from meth5 import MetH5File


def argtype_M5File(value):
    try:
        MetH5File(value, "r").get_chromosomes()
    except:
        raise argparse.ArgumentTypeError(f"Failed to read '{value}'. Is it a valid MetH5 file?")
    return Path(value)


def argtype_genomic_range(value: str):
    try:
        value = value.split(":")
        chrom = value[0]
        value = value[1].split("-")
        start = int(value[0])
        end = int(value[1])
        return dict(chrom=chrom, start=start, end=end)
    except:
        raise argparse.ArgumentTypeError(f"Failed to parse '{value}' as a genomic range. Must be chrom:start-end")
