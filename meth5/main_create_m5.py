from typing import List
import argparse
import tqdm
from pathlib import Path
from meth5 import MetH5File
from meth5.util import argtype_M5File

__description__ = "Create m5 file from Nanopolish result files"


def set_arguments(sc_args: argparse.ArgumentParser):
    sc_args.add_argument(
        "-i",
        "--input_paths",
        type=Path,
        nargs="+",
        required=True,
        help="Path(s) to Nanopolish result files or folder containing them",
    )
    
    sc_args.add_argument(
        "-o",
        "--output_file",
        type=Path,
        required=True,
        help="Output MetH5 file",
    )
    
    sc_args.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="No progress bar or warnings will be displayed",
    )
    
    sc_args.add_argument(
        "-c",
        "--compression",
        type=str,
        required=False,
        default="lzf",
        choices=["gzip", "lzf", "None"],
        help="Compression method for the MetH5 data structures. Use 'gzip' for smallest file size. Use 'lzf' for "
        "fastest access speed. For no compression you can also provide 'None', but is not recommended. "
        "Default: lzf",
    )
    
    sc_args.add_argument(
        "--allowed_chromosomes",
        type=str,
        nargs="+",
        required=False,
        default=None,
        help="Only include these chromosomes",
    )


def main(
    chunk_size: int,
    input_paths: List[Path],
    output_file: Path,
    compression: str,
    allowed_chromosomes: List[str],
    quiet: bool,
):
    if compression == "None":
        compression = None
    
    if allowed_chromosomes is not None and len(allowed_chromosomes) == 0:
        allowed_chromosomes = None
    
    input_files = []
    for input_path in input_paths:
        if not input_path.exists():
            raise ValueError(f"Cannot find path {input_path}")
        if input_path.is_file():
            input_files.append(input_path)
        elif input_path.is_dir():
            subfiles = list(input_path.iterdir())
            if len(subfiles) == 0:
                raise ValueError(f"Provided input path refers to a directory which is empty: {input_path}")
            input_files += [f for f in subfiles if f.is_file()]
    
    with MetH5File(output_file, chunk_size=chunk_size, mode="w", compression=compression) as m5_out:
        for input_file in tqdm.tqdm(input_files) if not quiet else input_files:
            m5_out.parse_and_add_nanopolish_file(
                input_file, postpone_sorting_until_close=True, include_chromosomes=allowed_chromosomes
            )
        m5_out.create_chunk_index()
