import argparse
from pathlib import Path
from typing import IO
import logging
from meth5.util import argtype_M5File

import pandas as pd
from meth5 import MetH5File

__description__ = "Annotate reads with read group"


def set_arguments(sc_args: argparse.ArgumentParser):
    sc_args.add_argument(
        "--m5file",
        required=True,
        type=argtype_M5File,
        help="MetH5 file containing methylation calls",
    )
    
    sc_args.add_argument(
        "--read_groups_key",
        type=str,
        required=True,
        help="Read groups key under which the groups should be stored",
    )
    
    sc_args.add_argument(
        "--read_group_file",
        type=Path,
        required=True,
        help="Tab-delimited file containing columns read_name and numeric group",
    )


def read_readgroups(readgroups_file: IO):
    """
    Reads file that assigns read to read groups (such as haplotypes, samples, clusters, etc)
    :param readgroups_file: path to the tab-separated file
    :return: pandas dataframe with columns "read_name" and "group"
    """
    # Loading
    try:
        read_groups = pd.read_csv(
            readgroups_file, sep="\t", header=0, dtype={"read_name": str, "group": str}, index_col=None
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
    all_groups = list(set(read_annotation["group"]))
    group_dict = {g:i for i,g in enumerate(all_groups)}
    read_annotation["group"] = read_annotation["group"].map(group_dict.get)
    group_dict = {v:k for k,v in group_dict.items()}
    read_annotation = read_annotation.set_index("read_name")["group"].to_dict()
    
    with MetH5File(m5file, mode="a", chunk_size=chunk_size) as m5:
        m5.annotate_read_groups(read_groups_key, read_annotation, labels=group_dict, exists_ok=True, overwrite=True)
