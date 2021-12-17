from typing import List
import tqdm
import pandas as pd
from pathlib import Path
from meth5.meth5 import MetH5File

def isint(s):
    try:
        _ = int(s)
        return True
    except:
        return False
    

def main(
    chunk_size: int,
    input_m5_files: List[Path],
    read_group_names: List[str],
    read_groups_key: str,
    output_file: Path,
    compression: str,
    allowed_chromosomes: List[str],
    quiet: bool,
):
    if compression == "None":
        compression = None
    
    if allowed_chromosomes is not None and len(allowed_chromosomes) == 0:
        allowed_chromosomes = None
    
    if len(read_group_names) != len(input_m5_files):
        raise ValueError(f"List of read group prefixes must match number of input files")
    
    with MetH5File(output_file, chunk_size=chunk_size, mode="w", compression=compression) as m5_out:

        rg_maps = {read_groups_key:{}}
        for i, input_file in enumerate(input_m5_files):
            with MetH5File(input_file, "r") as m5_in:
                print("Reading ", input_file)
                if "read_groups" in m5_in.h5_fp["reads"].keys():
                    for rg_key in m5_in.get_read_group_keys():
                        names = m5_in.get_all_read_groups(rg_key)
                        if rg_key not in rg_maps:
                            rg_maps[rg_key] = {}
                        read_group_ids = m5_in.h5_fp["reads"]["read_groups"][rg_key][()]
                        read_names = m5_in.h5_fp["reads"]["read_names_mapping"][()]
                        print("Reading read groups ", rg_key)
                        for rn, rg_id in zip(read_names,  tqdm.tqdm(read_group_ids)):
                            rg_maps[rg_key][rn.decode()] = names[rg_id]
                    
                    for rn in read_names:
                        rg_maps[read_groups_key][rn.decode()] = read_group_names[i]
                if allowed_chromosomes is None:
                    chromosomes = set(m5_in.get_chromosomes())
                else:
                    chromosomes = set(allowed_chromosomes).intersection(set(m5_in.get_chromosomes()))

                print("Copying methylation calls")
                with tqdm.tqdm(total=100, disable=quiet) as pbar:
                    progress = 0
                    percent_per_chrom = 100 / len(chromosomes)
                    for chromosome in chromosomes:
                        chrom_container = m5_in[chromosome]
                        percent_per_chunk = percent_per_chrom / chrom_container.get_number_of_chunks()
                        for chunk in chrom_container.get_chunk_ids():
                            values_container = chrom_container.get_chunk(chunk, overlap=False)
                            ranges = values_container.get_ranges()
                            llrs = values_container.get_llrs()
                            read_names = values_container.get_read_names()
                            df = pd.DataFrame(
                                {
                                    "chromosome": chromosome,
                                    "start": ranges[:, 0],
                                    "end": ranges[:, 1],
                                    "read_name": read_names,
                                    "log_lik_ratio": llrs,
                                }
                            )
                            m5_out.add_to_h5_file(df, postpone_sorting_until_close=True)
                            progress += percent_per_chunk
                            pbar.n = progress
                            pbar.refresh()
        print("Indexing")
        m5_out.create_chunk_index()
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
