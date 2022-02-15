from pathlib import Path
import argparse
from meth5 import MetH5File


def argtype_M5File(value):
    try:
        MetH5File(value, "r").get_chromosomes()
    except:
        raise argparse.ArgumentTypeError(f"Failed to read '{value}'. Is it a valid MetH5 file?")
    return Path(value)
