import time
from pathlib import Path
from typing import IO, Iterable, List
from multiprocessing import Queue, Process
import argparse
import logging

import tqdm
import pandas as pd
import numpy as np
from meth5.meth5 import MetH5File
from meth5.sparse_matrix import SparseMethylationMatrixContainer

from nanoepiseg.emissions import BernoulliPosterior
from nanoepiseg.hmm import SegmentationHMM
from nanoepiseg.postprocessing import cleanup_segmentation
from nanoepiseg.segments_csv_io import SegmentsWriterBED
from nanoepiseg.math import llr_to_p

def read_readgroups(readgroups_file: IO):
    """
    Reads file that assigns read to read groups (such as haplotypes, samples, clusters, etc)
    :param readgroups_file: path to the tab-separated file
    :return: pandas dataframe with columns "read_name" and "group"
    """
    # Loading
    try:
        read_groups = pd.read_csv(
            readgroups_file,
            sep="\t",
            header=0,
            dtype={"read_name": str, "group": int},
            index_col=None
        )
    except Exception as e:
        logging.error("Unable to read read groups file", e)
        raise e
    
    # Validation
    if len(read_groups.columns) == 2:
        should_colnames = ["read_name", "group"]
    else:
        raise ValueError(f"Invalid number of columns in read groups file (should be 2, was {len(read_groups.columns)})")
    
    if not all([col in read_groups.columns for col in should_colnames]):
        raise ValueError("Invalid column names in read groups file (should be %s)" % should_colnames.join(", "))
    
    return read_groups


def main(
    m5file: Path,
    read_groups_key: str,
    read_group_file: Path,
    chunk_size: int,
):
    read_annotation = read_readgroups(read_group_file)
    read_annotation = read_annotation.set_index("read_name")["group"].to_dict()
    
    with MetH5File(m5file, mode="a", chunk_size = chunk_size) as m5:
        m5.annotate_read_groups(read_groups_key, read_annotation, exists_ok = True, overwrite = True)
