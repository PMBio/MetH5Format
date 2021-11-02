from typing import List
import argparse
import tqdm
import pandas as pd
from pathlib import Path
from meth5.meth5 import MetH5File


def main(
    chunk_size: int,
    input_m5_files: List[Path],
    read_group_names: List[str],
    read_group_key: str,
    output_file: Path,
    compression: str,
    allowed_chromosomes: List[str],
    quiet: bool,
):
    if compression == "None":
        compression = None
    
    if allowed_chromosomes is not None and len(allowed_chromosomes) == 0:
        allowed_chromosomes = None
    
    if len(input_m5_files) == 0:
        raise ValueError(f"No input files found in input directory f{str(input_dir)}")
    
    if len(read_group_prefixes) != len(input_m5_files):
        raise ValueError(f"List of read group prefixes must match number of input files")
    
    with MetH5File(output_file, chunk_size=chunk_size, mode="w", compression=compression) as m5_out, tqdm.tqdm(
        total=100, disable=quiet
    ) as pbar:
        
        percent_per_file = 100.0 / len(input_m5_files)
        progress = 0
        rg_maps = {}
        for i, input_file in enumerate(input_m5_files):
            with MetH5File(input_file, "r") as m5_in:
                rg_names = {}
                for rg_key in m5_in.get_read_group_keys():
                    rg_names[rg_key] = m5_in.get_all_read_groups(rg_key)
                    print(rg_names)
                    if rg_key not in rg_maps:
                        rg_maps[rg_key] = {}
                
                if allowed_chromosomes is None:
                    chromosomes = set(m5_in.get_chromosomes())
                else:
                    chromosomes = set(allowed_chromosomes).intersection(set(m5_in.get_chromosomes()))
                
                percent_per_chrom = percent_per_file / len(chromosomes)
                for chromosome in chromosomes:
                    chrom_container = m5_in[chromosome]
                    percent_per_chunk = percent_per_chrom / chrom_container.get_number_of_chunks()
                    for chunk in chrom_container.get_chunk_ids():
                        values_container = chrom_container.get_chunk(chunk)
                        ranges = values_container.get_ranges()
                        read_names = values_container.get_read_names()
                        df = pd.DataFrame(
                            {
                                "chromosome": chromosome,
                                "start": ranges[:, 0],
                                "end": ranges[:, 1],
                                "read_name": read_names,
                                "log_lik_ratio": values_container.get_llrs(),
                            }
                        )
                        
                        m5_out.add_to_h5_file(df, postpone_sorting_until_close=True)
                        
                        for rg_key in rg_names.keys():
                            read_groups = values_container.get_read_groups(rg_key)
                            names = rg_names[rg_key]
                            for rn, rg_id in zip(read_names, read_groups):
                                rg_maps[rg_key][rn] = f"{read_group_prefixes[i]}_{names[rg_id]}"
                        
                        progress += percent_per_chunk
                        pbar.n = progress
                        pbar.refresh()
        
        m5_out.create_chunk_index()
        for rg_key in rg_maps:
            read_groups = list(set(rg_maps[rg_key].values()))
            labels_inv = {v: i for i, v in enumerate(read_groups)}
            rg_map = {rn: labels_inv[rg_label] for rn, rg_label in rg_maps[rg_key].items()}
            labels = {i: v for v, i in labels_inv.items()}
            m5_out.annotate_read_groups(rg_key, rg_map, labels)
