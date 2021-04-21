from __future__ import annotations
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import logging
from pathlib import Path
from typing import Union, List, Dict, IO, Tuple, Any
from types import FunctionType

from meth5.sparse_matrix import SparseMethylationMatrixContainer


def _unique_genomic_range(genomic_ranges: np.ndarray) -> np.ndarray:
    """Helper function to computes unique genomic ranges from a
    presorted list of ranges in linear time.

    :param genomic_ranges: Numpy array of shape (n, 2)
    :return: A subset of the input rows, containing only unique regions (m<n,2)
    """
    diff = np.ones_like(genomic_ranges)
    diff[1:] = genomic_ranges[1:] - genomic_ranges[:-1]
    idx = sorted(list(set(diff[:, 0].nonzero()[0]).intersection(set(diff[:, 1].nonzero()[0]))))
    return genomic_ranges[idx, :]


def compute_betascore(llrs, llr_threshold=2):
    num_llrs = (np.abs(llrs) > llr_threshold).sum()
    return (llrs > llr_threshold).sum() / num_llrs if num_llrs > 0 else np.nan

def create_sparse_matrix_from_samples(
    sample_met_containers: Dict[str, MethlyationValuesContainer], sample_prefix_readnames=False,
) -> SparseMethylationMatrixContainer:
    """Creates a SparseMethylationMatrixContainer from a dictionary of
    MethylationValuesContainer. Each key value pair represents one
    sample.

    This helper function can be used if data is stored as one Meth5 file per sample,
    but should be analyzed together.

    The resulting sparse matrix is stored as a csc_matrix and is  created directly to
    keep memory requirement low

    :param sample_met_containers: keys are sample names, values are the corresponding
    MethylationValuesContainers extracted from a MetH5File
    :param sample_prefix_readnames: If you are worried that read names between samples
    might contain duplicates (collisions), this will remove those by prefixing the
    read names with the sample name
    :return: SparseMethylationMatrixContainer where read_samples are set based on the
    dictionary keys of the input
    """
    # Decide on a canonical order of samples
    samples = list(sample_met_containers.keys())

    read_names_dict = {
        s: [
            s + r.decode() if sample_prefix_readnames else r.decode()
            for r in sample_met_containers[s].get_read_names_unique()
        ]
        for s in samples
    }
    genomic_ranges = {s: [r for r in sample_met_containers[s].get_ranges_unique()] for s in samples}

    # Creates a sample assignment for every methylation call
    sample_assignment = [s for s in samples for _ in read_names_dict[s]]
    read_names = [r for s in samples for r in read_names_dict[s]]
    genomic_ranges = [r for s in samples for r in genomic_ranges[s]]
    genomic_ranges = np.array(sorted(genomic_ranges, key=lambda x: x[0] * 10e10 + x[1]))
    genomic_ranges = _unique_genomic_range(genomic_ranges)

    coord_to_index_dict = {genomic_ranges[i, 0]: i for i in range(len(genomic_ranges))}

    read_dict = {read_names[i]: i for i in range(len(read_names))}

    # Fill out the actual data
    sparse_data = []
    sparse_x = []
    sparse_y = []
    for sample, llrs in sample_met_containers.items():
        range_ds = llrs.get_ranges()
        read_name_ds = llrs.get_read_names()

        llr_ds = llrs.get_llrs()
        sparse_data += list(llr_ds[:])
        sparse_x += [read_dict[sample + r.decode() if sample_prefix_readnames else r.decode()] for r in read_name_ds[:]]
        sparse_y += [coord_to_index_dict[gr[0]] for gr in range_ds[:]]

    # Create sparse matrix
    met_matrix = sp.csc_matrix((sparse_data, (sparse_x, sparse_y)))
    return SparseMethylationMatrixContainer(
        met_matrix, read_names, genomic_ranges[:, 0], genomic_ranges[:, 1], read_samples=sample_assignment,
    )


class MethlyationValuesContainer:
    """Manages access to the data (methylation llrs, read names, etc) of
    a genomic region on one chromosome."""

    def __init__(self, chromosome_container: ChromosomeContainer, start: int, end: int):
        """
        :param chromosome_container: Parent ChromosomeContainer object
        :param start: start index (not genomic location, but index in the dataframes)
        :param end: end index (not genomic location, but index in the dataframes)
        """
        self.chromosome = chromosome_container
        self.start = start
        self.end = end

    def get_read_names_unique(self) -> np.ndarray:
        """
        :return: Unique name of reads intersecting with this region
        """
        read_name_ds = self.chromosome.h5group["read_name"][self.start : self.end]
        reads_names, idx = np.unique(read_name_ds, return_index=True)
        reads_names = reads_names[np.argsort(idx)]
        return reads_names

    def get_ranges_unique(self) -> np.ndarray:
        """
        :return: Numpy array of shape (u, 2) for u unique genomic regions. Note that
        regions can overlap and can theoretically have the same starting but different
        end point. Ranges are sorted by start position first and then by end position.
        """
        ranges_ds = self.chromosome.h5group["range"][self.start : self.end]
        return _unique_genomic_range(ranges_ds)

    def get_ranges(self) -> np.ndarray:
        """
        :return: Numpy array of shape (n, 2) containing start and stop position of the
        genomic region of the associated methylation call.
        """
        return self.chromosome.h5group["range"][self.start : self.end, :]

    def get_llrs(self) -> np.ndarray:
        """
        :return: Numpy array of shape (n) containing the methylation call
        log-likelihood ratio
        """
        return self.chromosome.h5group["llr"][self.start : self.end]

    def get_read_names(self) -> np.ndarray:
        """
        :return: Numpy array of shape (n) containing the read name for each
        methylation call
        """
        return self.chromosome.h5group["read_name"][self.start : self.end]

    def get_read_groups(self, group_key: str) -> np.ndarray:
        """The Meth5 file can store multiple different groupings of
        methylation calls. Typically, this would be based on grouping
        reads (such as from read-phasing) but any sort of grouping of
        methylation calls is supported. This function returns the
        grouping of methylation calls.

        :param group_key: The group key under which the grouping has been stored
        :return: Numpy array of shape (n) containing the read group for each
        methylation call, given the grouping key.
        """
        return self.chromosome.h5group["read_groups"][group_key][self.start : self.end]

    def __compute_llr_site_aggregate(self, ranges, llrs, aggregation_fun):
        # Takes advantage of ranges being sorted
        range_diff = (np.diff(ranges[:, 0]) != 0) | (np.diff(ranges[:, 1]) != 0)
        # Changepoints where it goes from one range to the next
        range_cp = np.argwhere(range_diff).flatten() + 1
        range_start = [0, *range_cp]
        range_end = [*range_cp, llrs.shape[0]]
        # Calls aggregation function once for each unique range
        aggregated_llrs = np.array([aggregation_fun(llrs[rs:re]) for rs, re in zip(range_start, range_end)])
        return aggregated_llrs, ranges[range_start, :]

    def get_llr_site_aggregate(self, aggregation_fun: FunctionType) -> Tuple[np.ndarray, np.ndarray]:
        """Computes per-site an aggregate of the LLR. The provided
        aggregation function should take a numpy array and can return
        any arbitrary aggregate. The return value is a numpy array
        containing the aggregates for each unique genomic range.

        Note that ranges with same same startpoint but different endpoint will
        be considered as two separate ranges

        :param aggregation_fun: Function that takes a numpy array of llrs and returns the aggregate

        :return: Tuple consisting of:
          * aggregated llrs
          * ranges for each aggregation
        """
        llrs = self.get_llrs()
        ranges = self.get_ranges()

        return self.__compute_llr_site_aggregate(ranges, llrs, aggregation_fun)

    def get_llr_site_median(self):
        """Calls get_llr_site_aggregate with np.median as an aggregation function"""
        return self.get_llr_site_aggregate(np.median)

    def get_llr_site_rate(self, llr_threshold=2):
        """Calls get_llr_site_aggregate computing the methylation betascore"""
        return self.get_llr_site_aggregate(lambda llrs: compute_betascore(llrs, llr_threshold))
    
    def get_llr_read_aggregate(self, aggregation_fun: FunctionType) -> Dict[str, Any]:
        """Computes per-read an aggregate of the LLR. The provided
        aggregation function should take a numpy array and can return
        any arbitrary aggregate. The return value is a numpy array
        containing the aggregates for each unique genomic range.

        Note that ranges with same same startpoint but different endpoint will
        be considered as two separate ranges

        :param aggregation_fun: Function that takes a numpy array of llrs and returns the aggregate

        :return: Tuple consisting of:
          * aggregation result
          * ranges for each aggregation
        """
        llrs = self.get_llrs()
        reads = self.get_read_names()
        readset = set(reads)

        aggregated_llrs = {read: aggregation_fun(llrs[reads == read]) for read in readset}
        return aggregated_llrs


    def get_llr_site_readgroup_aggregate(
        self, group_key: str, aggregation_fun: FunctionType
    ) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """For each read group, computes a per-site aggregate of the LLR. The provided
        aggregation function should take a numpy array and can return
        any arbitrary aggregate. The return value is a dictionary with the key
        being each read group and the value being a tuple with the numpy arrays
        containing the aggregates for each range and in second position the genomic ranges

        Note that ranges with same same startpoint but different endpoint will
        be considered as two separate ranges

        :param group_key: The group key under which the grouping has been stored
        :param aggregation_fun: Function that takes a numpy array of llrs and returns the aggregate

        :return: {readgroup_key: (aggregated llrs, ranges for each aggregation)
        """
        all_llrs = self.get_llrs()
        all_ranges = self.get_ranges()
        all_groups = self.get_read_groups(group_key)

        return {
            group: self.__compute_llr_site_aggregate(
                all_ranges[all_groups == group], all_llrs[all_groups == group], aggregation_fun
            )
            for group in set(all_groups)
        }

    def get_llr_site_readgroup_rate(self, group_key: str, llr_threshold: float = 2):
        """Calls get_llr_site_readgroup_aggregate computing the methylation betascore"""
        return self.get_llr_site_readgroup_aggregate(group_key, lambda llrs: compute_betascore(llrs, llr_threshold))

    def to_sparse_methylation_matrix(self, read_groups_key: str = None) -> SparseMethylationMatrixContainer:
        """Creates a SparseMethylationMatrixContainer from the values in
        this container. If a read_groups_key is provided, then Meth5
        file will be checked for a matching read group annotation, which
        will then serve to define the samples in the
        SparseMethylationMatrixContainer.

        The resulting sparse matrix is stored as a csc_matrix and is  created
        directly to keep memory requirement low

        :param read_groups_key: The key in the Meth5 file under which the read groups
        (samples) can be found
        :return: SparseMethylationMatrixContainer or None
        """
        # Define canonical order of read names
        read_names = [r.decode() for r in self.get_read_names_unique()]
        genomic_ranges = self.get_ranges_unique()

        # Assigns y coordinate in the matrix to a genomic position
        coord_to_index_dict = {genomic_ranges[i, 0]: i for i in range(len(genomic_ranges))}

        # Assigns x coordinate in the matrix to a read name
        read_dict = {read_names[i]: i for i in range(len(read_names))}

        read_name_ds = self.get_read_names()

        sparse_data = self.get_llrs()[:]
        sparse_x = [read_dict[r.decode()] for r in read_name_ds]
        sparse_y = [coord_to_index_dict[p] for p in self.get_ranges()[:, 0]]

        if read_groups_key is not None:
            read_groups_ds = self.get_read_groups(read_groups_key)
            read_samples_dict = {rn.decode(): rg for (rn, rg) in zip(read_name_ds[:], read_groups_ds[:])}
            read_samples = np.array([read_samples_dict[r] for r in read_names])
        else:
            read_samples = None

        """Note: It's important to provide "shape" in the constructor, in case
        the matrix is empty. Otherwise the csc_matrix constructor will raise
        an error for not being able to infer the dimensions of the matrix"""
        met_matrix = sp.csc_matrix((sparse_data, (sparse_x, sparse_y)), shape=(len(read_names), len(genomic_ranges)))
        return SparseMethylationMatrixContainer(
            met_matrix, read_names, genomic_ranges[:, 0], genomic_ranges[:, 1], read_samples=read_samples,
        )


class ChromosomeContainer:
    """Manages access to the data of a single chromosome and provides
    functions for efficient subsetting (e.g. by chunk or by genomic
    region)"""

    def __init__(self, chromosome_group: h5py.Group, chunk_size: int):
        """
        :param chromosome_group: h5py.Group object inside the Meth5 file containing
        values for this chromosome
        :param chunk_size: chunk size to use for hdf5 dataframes
        """
        self.h5group = chromosome_group
        self.chunk_size = chunk_size

    def __len__(self) -> int:
        """
        :return: number of methylation calls on the entire chromosome
        """
        return len(self.h5group["range"])

    def get_number_of_chunks(self) -> int:
        """
        :return: given length and chunk size, returns the number of chunks
        """
        num_chunks = len(self) // self.chunk_size
        if len(self) % self.chunk_size != 0:
            num_chunks += 1
        return num_chunks

    def get_chunk_ids(self) -> List[int]:
        """
        :return: List of integer ids, one for each chunk.
        In the current implementation it's just a running counter
        """
        return [i for i in range(self.get_number_of_chunks())]

    def _seek_overlap_ranges_backwards(self, chunk_id: int, start_value: int = -1) -> int:
        """This helper function recursively looks backwards starting
        from a specified chunk, and returns the index of the first
        position in the dataframes that contains a methylation call for
        the same genomic site as the start of the provided chunk. Used
        to make sure all methylation calls (from all reads) are
        included.

        :param chunk_id: starting chunk id
        :param start_value: used in recursion only - don't overwrite it
        :return: first index for included sites
        """
        last = min(len(self), self.chunk_size * (chunk_id + 1)) - 1
        if start_value == -1:
            start_value = self.h5group["range"][self.chunk_size * chunk_id, 0]

        starts = self.h5group["range"][(self.chunk_size * chunk_id) : last, 0]
        matches = np.arange(len(starts))[starts == start_value]

        if len(matches) == 0:
            # Nothing in this chunk, return beginning of the chunk we came from
            return self.chunk_size * (chunk_id + 1)

        if matches[0] == 0 and chunk_id > 0:
            # All of this chunk is the same range, we need to go deeper
            return self._seek_overlap_ranges_backwards(chunk_id - 1, start_value=start_value)

        # Part of this chunk has entries for this start position
        return self.chunk_size * chunk_id + matches[0]

    def _seek_overlap_ranges_forwards(self, chunk_id, end_value=-1):
        """This helper function recursively looks forwards starting from
        the end of a specified chunk, and returns the index of the last
        position in the dataframes that contains a methylation call for
        the same genomic site as the end of the provided chunk. Used to
        make sure all methylation calls (from all reads) are included.

        :param chunk_id: starting chunk id
        :param end_value: used in recursion only - don't overwrite it
        :return: last index for included sites
        """
        last = min(len(self), self.chunk_size * (chunk_id + 1)) - 1

        if end_value == -1:
            end_value = self.h5group["range"][last, 0]

        ends = self.h5group["range"][(self.chunk_size * chunk_id) : (last + 1), 0]

        matches = np.arange(len(ends))[ends == end_value]

        if len(matches) == 0:
            # Nothing in this chunk, return end of the chunk we came from
            return self.chunk_size * chunk_id - 1

        if matches[-1] == self.chunk_size - 1 and chunk_id < self.get_number_of_chunks() - 1:
            # All of this chunk is the same range, we need to go deeper
            return self._seek_overlap_ranges_forwards(chunk_id + 1, end_value=end_value)

        # Part of this chunk has entries for this end position
        return self.chunk_size * chunk_id + matches[-1]

    def get_chunk(self, chunk_id: int, overlap=True) -> MethlyationValuesContainer:
        """Returns a MethlyationValuesContainer providing access to the
        values of said chunk, and, if overlap=True, includes values of
        neighboring chunks if they are in the same genomic ranges, such
        as to avoid having a subset of reads of one location in one
        chunk and the rest in the other.

        :param chunk_id: The chunk id (see get_chunk_ids)
        :param overlap: Whether to look for same-region locations in
        neighboring chunks
        :return: MethlyationValuesContainer
        """
        if overlap:
            earliest_pos = self._seek_overlap_ranges_backwards(chunk_id)
            latest_pos = self._seek_overlap_ranges_forwards(chunk_id) + 1
        else:
            earliest_pos = self.chunk_size * chunk_id
            latest_pos = min(self.chunk_size * (chunk_id + 1), len(self))

        return MethlyationValuesContainer(self, earliest_pos, latest_pos)

    def create_chunk_index(self, force_update=False):
        """Needs Meth5 file to be open in write or append mode. Creates
        an additional datastructure in the HDF5 file that stores an
        index that stores genomic start and end site of a chunk, for
        fast searching.

        :param force_update: Whether an existing index should be overwritten
        """

        if "chunk_ranges" in self.h5group.keys() and not force_update:
            return

        index = np.zeros((self.get_number_of_chunks(), 2))
        num_ranges = self.h5group["range"].shape[0]
        for chunk_i, start_i in enumerate(range(0, num_ranges, self.chunk_size)):
            end_i = min(num_ranges - 1, start_i + self.chunk_size)
            index[chunk_i, 0] = self.h5group["range"][start_i, 0]
            index[chunk_i, 1] = self.h5group["range"][end_i, 1]

        if "chunk_ranges" in self.h5group.keys():
            self.h5group["chunk_ranges"].resize(index.shape)
            self.h5group["chunk_ranges"][:] = index
        else:
            self.h5group.create_dataset(name="chunk_ranges", data=index, dtype=int, maxshape=(None, 2))
        self.h5group.attrs["chunk_size"] = self.chunk_size

    def get_all_values(self) -> MethlyationValuesContainer:
        """Returns a MethlyationValuesContainer providing access to all sites on the chromosome
        Very inefficient and therefore not recommended. Chunk-based operations are recommended.
        :return: MethlyationValuesContainer
        """
        return MethlyationValuesContainer(self, 0, self.h5group["range"].shape[0])

    def get_values_in_range(self, genomic_start: int, genomic_end: int) -> MethlyationValuesContainer:
        """Returns a MethlyationValuesContainer providing access to the
        specified genomic region.

        Needs an index created by create_chunk_index, since the chunk
        index is used for fast searching.
        :param genomic_start: Genomic start location on the chromosome
        :param genomic_end: Genomic end location on the chromosome
        :return: MethlyationValuesContainer or None if no values for the region are
         available
        """
        if "chunk_size" not in self.h5group.attrs.keys():
            raise ValueError("Random access to ranges only allowed if index exists. Call create_chunk_index")
        index_chunk_size = self.h5group.attrs["chunk_size"]
        index = self.h5group["chunk_ranges"][:]

        # First find the right chunk for start and end
        chunk_indices = np.arange(len(index))[(index[:, 0] < genomic_end) & (genomic_start <= index[:, 1])]

        if len(chunk_indices) == 0:
            # If no chunk contains these values
            return MethlyationValuesContainer(self, 0, 0)

        start_chunk = chunk_indices[0]
        end_chunk = chunk_indices[-1]

        # Find precise start point
        start_index = start_chunk * index_chunk_size
        start_chunk_start = start_chunk * index_chunk_size
        start_chunk_end = min(len(self) - 1, (start_chunk + 1) * index_chunk_size)
        start_chunk_ranges = self.h5group["range"][start_chunk_start:start_chunk_end, :]
        start_in_range_indices = np.arange(len(start_chunk_ranges))[start_chunk_ranges[:, 1] >= genomic_start]
        if len(start_in_range_indices) > 0:
            # Add index of first value that is in the range
            start_index += start_in_range_indices[0]

        # Find precise end point
        end_index = end_chunk * index_chunk_size
        end_chunk_start = end_chunk * index_chunk_size
        end_chunk_end = min(len(self) - 1, (end_chunk + 1) * index_chunk_size)
        end_chunk_ranges = self.h5group["range"][end_chunk_start:end_chunk_end, :]
        end_oor_indices = np.arange(len(end_chunk_ranges))[end_chunk_ranges[:, 0] >= genomic_end]
        if len(end_oor_indices) > 0:
            # Add index of first value that is out of range
            end_index += end_oor_indices[0]
        else:
            # If all values in the chunk are in the range
            end_index = min(len(self), end_index + index_chunk_size)

        return MethlyationValuesContainer(self, start_index, end_index)


class MetH5File:
    """Main wrapper for Meth5 files."""

    def __init__(self, h5filepath: Union[str, Path, IO], mode: str = "r", chunk_size=int(10e5), compression="gzip"):
        """Initializes Meth5File and directly opens the file pointer.

        :param h5filepath: Path to Meth5 file or IO object providing access to it
        :param mode: h5py.File mode (typically "r", "w", or "a")
        :param chunk_size: chunk size to be used for HDF5 dataframes as well as for
        indexing and searching
        """
        self.h5filepath = h5filepath
        self.mode = mode
        self.chunk_size = chunk_size
        self.h5_fp: h5py.File = None
        self.chrom_container_cache = {}
        self.log = logging.getLogger("NET:MetH5")
        self.compression = compression
        self.h5_fp = h5py.File(self.h5filepath, mode=self.mode)

    def __enter__(self):
        return self

    def close(self):
        """Close HDF file pointer."""
        self.h5_fp.close()

    def resort_chromosome(self, chrom: str):
        """Forces resorting values of one chromosome by range"""
        chrom_group = self.h5_fp["chromosomes"][chrom]
        sort_order = np.argsort(chrom_group["range"][:, 0], kind="mergesort")
        logging.debug("Re-sorting h5 entries for chromosome %s" % chrom)
        chrom_group["range"][:] = np.array(chrom_group["range"])[sort_order]
        chrom_group["llr"][:] = np.array(chrom_group["llr"])[sort_order]
        chrom_group["read_name"][:] = np.array(chrom_group["read_name"])[sort_order]
        chrom_group.attrs["is_sorted"] = True

    def resort_unsorted_chromosomes(self):
        """Resorts only those chromosomes that are unsorted (have the "is_sorted" attribute set to False)"""
        for chrom in self.get_chromosomes():
            if not self.h5_fp["chromosomes"][chrom].attrs.get("is_sorted", True):
                self.resort_chromosome(chrom)

    def __exit__(self, exittype, exitvalue, traceback):
        try:
            self.resort_unsorted_chromosomes()
        except:
            pass
        self.close()

    def _create_or_extend(self, parent_group: h5py.Group, name: str, shape: Tuple, data: np.ndarray, **kwargs):
        """Internal helper function that either creates a dataframe if
        it doesn't exist or it extends it by using the h5py resize
        function.

        :param parent_group: chromosome group
        :param name: name of the dataframe in the group
        :param shape: shape of the added data (not new shape after extending)
        :param data: data to be added
        :param kwargs: passed on to create_dataset only if it doesn't exist
        """
        if name not in parent_group.keys():
            parent_group.create_dataset(name=name, shape=shape, data=data, **kwargs)
        else:
            num_data = data.shape[0] if hasattr(data, "shape") else len(data)
            ds = parent_group[name]
            old_shape = ds.shape
            new_shape = (old_shape[i] + (num_data if i == 0 else 0) for i in range(len(old_shape)))
            ds.resize(new_shape)

            ds[old_shape[0] :] = data

            self.log.debug("Extended from %s to %s" % (old_shape, ds.shape))

    def add_to_h5_file(
        self, cur_df: pd.DataFrame, include_chromosomes: List[str] = None, postpone_sorting_until_close=False
    ):
        """Add data from a pandas Dataframe which is the result of
        reading a nanopolish output file. Must at least contain the
        columns "chromosome", "read_name", "start", "end",
        "log_lik_ratio".

        :param cur_df: pandas dataframe containing nanopolish output
        :param include_chromosomes: List of chromosome names to be included. Recommended
        if your mapping contains lots of alternative contigs that you don't plan to use
        downstream anyways. Can greatly improve performance. If None, all chromosomes are included.
        Default: None
        """
        main_group = self.h5_fp.require_group("chromosomes")

        for chrom in set(cur_df["chromosome"].astype("str")):
            if include_chromosomes is not None and chrom not in include_chromosomes:
                continue
            self.log.debug("Adding sites from chromosome %s to h5 file" % chrom)

            chrom_calls = cur_df.loc[cur_df["chromosome"] == chrom]
            print(f"Adding {chrom_calls.shape[0]} from chrom {chrom}")
            n = chrom_calls.shape[0]
            read_names = [read.encode() for read in chrom_calls["read_name"]]
            read_name_len = len(read_names[0])
            assert all([len(read) for read in read_names])

            chrom_group = main_group.require_group(chrom)

            self._create_or_extend(
                parent_group=chrom_group,
                name="range",
                shape=(n, 2),
                dtype=int,
                data=chrom_calls[["start", "end"]],
                compression=self.compression,
                chunks=(self.chunk_size, 2),
                maxshape=(None, 2),
            )
            # TODO Add strand as a (bool) dataframe
            self._create_or_extend(
                parent_group=chrom_group,
                name="llr",
                shape=(n,),
                dtype=float,
                data=chrom_calls["log_lik_ratio"],
                compression=self.compression,
                chunks=(self.chunk_size,),
                maxshape=(None,),
            )
            # TODO Switch to indexing reads numerically and storing a map of read_names
            # to index in another dataframe (or attr)
            self._create_or_extend(
                parent_group=chrom_group,
                name="read_name",
                shape=(n,),
                dtype="S%d" % read_name_len,
                data=read_names,
                compression=self.compression,
                chunks=(self.chunk_size,),
                maxshape=(None,),
            )

            # TODO think of a way to do this that doesn't require loading one entire
            # dataset into memory
            if postpone_sorting_until_close:
                chrom_group.attrs["is_sorted"] = False
            else:
                self.resort_chromosome(chrom)

    def parse_and_add_nanopolish_file(self, nanopolish_file: Union[str, Path], read_chunk_lines=1e6, **kwargs):
        """Reads nanopolish output file and appends data to the Meth5
        file.

        :param nanopolish_file: Path to nanopolish file
        :param include_chromosomes: List of chromosome names to be included. Recommended
        if your mapping contains lots of alternative contigs that you don't plan to use
        downstream anyways. Can greatly improve performance. If None, all chromosomes are included.
        Default: None
        """
        cur_df = pd.read_csv(nanopolish_file, sep="\t", dtype={"chromosome": str}, chunksize=read_chunk_lines)
        for df_chunk in cur_df:
            self.add_to_h5_file(df_chunk, **kwargs)

    def get_chromosomes(self) -> List[str]:
        """
        :return: unsorted list of chromosomes in Meth5File
        """

        return [str(k) for k in self.h5_fp["chromosomes"].keys()]

    def __getitem__(self, chromosome: str) -> ChromosomeContainer:
        """Returns ChromosomeContainer object managing access to values
        of the given chromosome.

        :param chromosome: the chromosome name
        :return: ChromosomeContainer object
        """

        if chromosome not in self.h5_fp["chromosomes"].keys():
            return None

        if not self.h5_fp["chromosomes"][chromosome].attrs.get("is_sorted", True):
            raise ValueError(
                "MetH5 file has been manipulated and sorting has been postponed. Need to resort before"
                "accessing values."
            )

        if chromosome in self.chrom_container_cache.keys():
            return self.chrom_container_cache[chromosome]
        else:
            ret = ChromosomeContainer(self.h5_fp["chromosomes"][chromosome], self.chunk_size)
            self.chrom_container_cache[chromosome] = ret
            return ret

    def create_chunk_index(self, *args, **kwargs):
        """Create chunk index for each chromosome. Also performs resorting if necessary.

        See documentation of ChromosomeContainer.create_chunk_index
        """
        self.resort_unsorted_chromosomes()
        for chromosome in self.get_chromosomes():
            self[chromosome].create_chunk_index(*args, **kwargs)

    # TODO Implement a method that returns methylation values for one read group

    def annotate_read_groups(
        self, read_group_key: str, map: Dict[str, int], labels: Dict[int, str] = None, exists_ok=False, overwrite=False
    ):
        """Store read group annotation in the Meth5 file, which can
        later be accessed through the MethylationValuesContainer object.

        Since Meth5 format can store multiple read group annotations, a key for each
        annotion is needed.
        :param read_group_key: key under which this annotation should be stored
        :param map: maps read names to read group
        :param labels: maps the read group key (int) to a readable label (string)
        :param exists_ok: if False, an error will be thrown if a grouping with this key
        already exists (default=False)
        :param overwrite: if exists_ok=True and overwrite=True, an existing mapping will
        be updated. If exists_ok=True and overwrite=False, nothing will be done in case
        a grouping with this key already exists
        """
        for chromosome in self.get_chromosomes():
            chr_g = self.h5_fp["chromosomes"][chromosome]
            rg_g = chr_g.require_group("read_groups")
            if read_group_key in rg_g.keys():
                if not exists_ok:
                    raise ValueError("Cannot annotate read groups - group assignment with this key " "already exists")
                elif not overwrite:
                    continue

            rg_assignment = [map.get(read.decode(), -1) for read in chr_g["read_name"][:]]
            rg_ds = rg_g.require_dataset(name=read_group_key, dtype=int, shape=(len(rg_assignment),), maxshape=(None,),)
            rg_ds[:] = rg_assignment

            rg_ds.attrs.clear()
            if labels is not None:
                # TODO: Find a nicer way to do this:
                # I originally intended to store the int keys, but hdf5 doesnt support
                # integers as keys in attributes dictionary....
                labels = {str(k): v for k, v in labels.items()}
                rg_ds.attrs.update(labels)
