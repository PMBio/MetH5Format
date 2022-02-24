from typing import List, Dict
import argparse
from typing import List
import tqdm
import pandas as pd
from pathlib import Path
from meth5 import MetH5File
from meth5.util import argtype_M5File


def isint(s):
    try:
        _ = int(s)
        return True
    except:
        return False


__description__ = "Merge m5 file from Nanopolish result files"


def set_arguments(sc_args: argparse.ArgumentParser):
    sc_args.add_argument(
        "--input_m5_files",
        type=argtype_M5File,
        required=True,
        nargs="+",
        help="List of MetH5 files",
    )
    
    group_rg = sc_args.add_mutually_exclusive_group(required=True)
    group_rg.add_argument(
        "--read_group_names",
        type=str,
        nargs="+",
        help="One name per input file",
    )
    group_rg.add_argument(
        "--no_read_groups", action="store_true", help="No read groups in output file",
    )

    sc_args.add_argument("--no_read_names", action="store_true", help="No read names in output file", )

    sc_args.add_argument(
        "--read_groups_key",
        type=str,
        required=True,
        help="Read groups key under which the groups should be stored",
    )
    
    sc_args.add_argument(
        "--output_file",
        type=Path,
        required=True,
        help="Output MetH5 file",
    )
    
    sc_args.add_argument(
        "--quiet",
        action="store_true",
        help="No progress bar or warnings will be displayed",
    )
    
    sc_args.add_argument(
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


def compute_total_chrom_sizes(input_m5_files: List[str]) -> Dict[str, int]:
    chrom_size = {}
    for i, input_file in enumerate(input_m5_files):
        with MetH5File(input_file, "r") as m5_in:
            for chrom in m5_in.get_chromosomes():
                if chrom not in chrom_size:
                    chrom_size[chrom] = 0
                chrom_size[chrom] += len(m5_in[chrom])
    return chrom_size


def main(
    chunk_size: int,
    input_m5_files: List[Path],
    read_group_names: List[str],
    read_groups_key: str,
    output_file: Path,
    compression: str,
    allowed_chromosomes: List[str],
    no_read_groups: bool,
    quiet: bool,
    no_read_names: bool,
):
    if compression == "None":
        compression = None
    
    if allowed_chromosomes is not None and len(allowed_chromosomes) == 0:
        allowed_chromosomes = None
    
    if read_groups_key is not None:
        if len(read_group_names) != len(input_m5_files):
            raise ValueError(f"List of read group prefixes must match number of input files")
    
    total_chrom_sizes = compute_total_chrom_sizes(input_m5_files)
    
    with MetH5File(
        output_file,
        chunk_size=chunk_size,
        mode="w",
        compression=compression,
        max_calls=total_chrom_sizes,
    ) as m5_out:
        
        if read_groups_key is not None:
            rg_maps = {read_groups_key: {}}
        else:
            rg_maps = {}
        
        read_id_offset = 0
        for i, input_file in enumerate(input_m5_files):
            max_read_id_local = 0
            with MetH5File(input_file, "r") as m5_in, tqdm.tqdm(total=100, disable=quiet) as pbar:
                print("Reading ", input_file)
                read_read_groups(m5_in, no_read_groups, read_group_names[i], read_groups_key, rg_maps)
                
                if allowed_chromosomes is None:
                    chromosomes = set(m5_in.get_chromosomes())
                else:
                    chromosomes = set(allowed_chromosomes).intersection(set(m5_in.get_chromosomes()))
                
                print("Copying methylation calls")
                
                progress = 0
                percent_per_chrom = 100 / len(chromosomes)
                for chromosome in chromosomes:
                    chrom_container = m5_in[chromosome]
                    percent_per_chunk = percent_per_chrom / chrom_container.get_number_of_chunks()
                    for chunk in chrom_container.get_chunk_ids():
                        values_container = chrom_container.get_chunk(chunk, overlap=False)
                        ranges = values_container.get_ranges()
                        llrs = values_container.get_llrs()
                        
                        if no_read_names:
                            read_names = values_container.get_read_ids()
                            read_names += read_id_offset
                            max_read_id_local = max(max_read_id_local, max(read_names))
                            read_names_key = "read_id"
                        else:
                            read_names = values_container.get_read_names()
                            read_names_key = "read_name"
                        
                        df = pd.DataFrame(
                            {
                                "chromosome": chromosome,
                                "start": ranges[:, 0],
                                "end": ranges[:, 1],
                                read_names_key: read_names,
                                "log_lik_ratio": llrs,
                            }
                        )
                        m5_out.add_to_h5_file(df, postpone_sorting_until_close=True)
                        progress += percent_per_chunk
                        pbar.n = progress
                        pbar.refresh()
        print("Indexing")
        m5_out.create_chunk_index()
        write_read_groups(m5_out, no_read_groups, rg_maps)


def write_read_groups(m5_out, no_read_groups, rg_maps):
    if not no_read_groups:
        print("Writing read groups")
        for rg_key in rg_maps:
            read_groups = list(set(rg_maps[rg_key].values()))
            if all(isint(x) for x in read_groups):
                labels_inv = {v: int(v) for v in read_groups}
            else:
                labels_inv = {v: i for i, v in enumerate(read_groups)}
            
            rg_map = {rn: labels_inv[rg_label] for rn, rg_label in rg_maps[rg_key].items()}
            labels = {i: v for v, i in labels_inv.items()}
            m5_out.annotate_read_groups(rg_key, rg_map, labels)


def read_read_groups(m5_in, no_read_groups, read_group_name, read_groups_key, rg_maps):
    if "read_groups" in m5_in.h5_fp["reads"].keys() and not no_read_groups:
        for rg_key in m5_in.get_read_group_keys():
            names = m5_in.get_all_read_groups(rg_key)
            if rg_key not in rg_maps:
                rg_maps[rg_key] = {}
            read_group_ids = m5_in.h5_fp["reads"]["read_groups"][rg_key][()]
            read_names = m5_in.h5_fp["reads"]["read_names_mapping"][()]
            print("Reading read groups ", rg_key)
            for rn, rg_id in zip(read_names, tqdm.tqdm(read_group_ids)):
                rg_maps[rg_key][rn.decode()] = names.get(rg_id, "unknown")
        
        if read_groups_key is not None:
            for rn in read_names:
                rg_maps[read_groups_key][rn.decode()] = read_group_name
