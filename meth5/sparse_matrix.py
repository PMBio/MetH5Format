from __future__ import annotations
from typing import List, Tuple, Union
import numpy as np
import scipy.sparse as sp


class SparseMethylationMatrixContainer:
    """This class manages a sparse matrix with dimensions
    (number_of_reads, number_of_sites).

    It is designed to hold log-likelihood ratios (LLR). Since an LLR of
    0 means a 50% chance of methylation (i.e. no information), a sparse
    matrix is particularly well suited to store them. Missing
    methylation calls are thus simply treated the same as a VERY
    uncertain methylation call.
    """
    
    def __init__(
        self,
        met_matrix: sp.csc_matrix,
        read_names: np.ndarray,
        genomic_coord_start: np.ndarray,
        genomic_coord_end: np.ndarray,
        read_samples: np.ndarray = None,
    ):
        """
        :param met_matrix: scipy csc_matrix of shape (n_reads, n_sites)
        :param read_names: numpy array of shapoe (n_reads) containing read names
        for each x-axis index
        :param genomic_coord_start: numpy array of shape (n_sites) containing
        genomic start coordinates for each y-axis index
        :param genomic_coord_end: numpy array of shape (n_sites) containing g
        enomic end coordinates for each y-axis index (e.g. when calls are for
        subsequences instead of basepair resolution)
        :param read_names: numpy array of shapoe (n_reads) containing read sample
        assignment for each x-axis index (can be None if no read grouping is required)
        """
        if met_matrix.shape[0] != len(read_names):
            raise ValueError(
                "List of read names must match dimension 0 of matrix (%d != %d)"
                % (met_matrix.shape[0], len(read_names))
            )
        if met_matrix.shape[1] != len(genomic_coord_start):
            raise ValueError(
                "Genomic start coordinates must match dimension 1 of matrix (%d != %d)"
                % (met_matrix.shape[1], len(genomic_coord_start))
            )
        if met_matrix.shape[1] != len(genomic_coord_end):
            raise ValueError(
                "Genomic end coordinates must match dimension 1 of matrix (%d != %d)"
                % (met_matrix.shape[1], len(genomic_coord_end))
            )
        if read_samples is not None:
            if met_matrix.shape[0] != len(read_samples):
                raise ValueError(
                    "Read sample assignment must match dimension 0 of matrix (%d != %d)"
                    % (met_matrix.shape[0], len(read_samples))
                )
        
        self.met_matrix = met_matrix
        self.shape = self.met_matrix.shape
        self.read_names = np.array(read_names)
        self.genomic_coord = np.array(genomic_coord_start)
        self.genomic_coord_end = np.array(genomic_coord_end)
        self.coord_to_index_dict = {genomic_coord_start[i]: i for i in range(len(genomic_coord_start))}
        self.read_samples = read_samples
    
    def _compact(self):
        """Helper function to be called after the matrix has been
        subset.

        Iteratively removes genomic sites for which there are no
        methylation calls, and then reads for which there are no calls,
        until all reads and sites have at least one non-zero methylation
        call
        """
        while True:
            read_has_values = np.array(((self.met_matrix != 0).sum(axis=1) > 0)).flatten()
            site_has_values = np.array(((self.met_matrix != 0).sum(axis=0) > 0)).flatten()
            self.met_matrix = self.met_matrix[read_has_values, :][:, site_has_values]
            self.read_names = self.read_names[read_has_values]
            self.genomic_coord = self.genomic_coord[site_has_values]
            self.genomic_coord_end = self.genomic_coord_end[site_has_values]
            if self.read_samples is not None:
                self.read_samples = self.read_samples[read_has_values]
            
            if (~read_has_values).sum() == 0 and (~site_has_values).sum() == 0:
                break
        self.shape = self.met_matrix.shape
        self.coord_to_index_dict = {self.genomic_coord[i]: i for i in range(len(self.genomic_coord))}
    
    def get_submatrix(self, start: int, end: int, compact: bool = True) -> SparseMethylationMatrixContainer:
        """Creates a submatrix containing only the genomic positions
        between start and end index.

        Note: start and end are matrix indices, not genomic coordinates
        Use get_submatrix_from_genomic_locations if you want to instead subset
        the matrix based on genomic coordinates.

        :param start: start index in matrix
        :param end: end index in matrix
        :param compact: whether or not to reads that no longer have any data
        :return: new SparseMethylationMatrixContainer object
        """
        sub_met_matrix = self.met_matrix[:, start:end]
        sub_genomic_coord = self.genomic_coord[start:end]
        sub_genomic_coord_end = self.genomic_coord_end[start:end]
        
        ret = SparseMethylationMatrixContainer(
            sub_met_matrix,
            self.read_names,
            sub_genomic_coord,
            sub_genomic_coord_end,
            read_samples=self.read_samples,
        )
        if compact:
            ret._compact()
        return ret
    
    def get_submatrix_from_genomic_locations(
        self, start_base: int, end_base: int, compact: bool = True
    ) -> SparseMethylationMatrixContainer:
        """Returns a new SparseMethylationMatrixContainer containing
        only data for the provided genomic reigon.

        :param start_base: genomic start location
        :param end_base: genomic end location
        :return: new SparseMethylationMatrixContainer object
        """
        start = self.coord_to_index_dict[start_base]
        end = self.coord_to_index_dict[end_base]
        return self.get_submatrix(start, end, compact=compact)
    
    def get_submatrix_from_genomic_locations_mask(
        self, allowed_sites: np.ndarray(), compact: bool = True
    ) -> SparseMethylationMatrixContainer:
        """Creates a submatrix containing only the given sites and their
        methylation calls.

        Compatibility note: If used with scipy 1.5.2, allowed_sites can be
        any iterable. If used with scipy 1.4.2, allowed_sites must be a numpy
        array. If providing a list, it may sometimes throw ValueErrors

        :param allowed_sites: boolean array of shape (n_sites)
        :param compact: whether or not to genomic sites that no longer have
        any data
        :return: new SparseMethylationMatrixContainer object
        """
        
        sub_met_matrix = self.met_matrix[:, allowed_sites]
        sub_genomic_coord = self.genomic_coord[allowed_sites]
        sub_genomic_coord_end = self.genomic_coord_end[allowed_sites]
        
        ret = SparseMethylationMatrixContainer(
            sub_met_matrix,
            self.read_names,
            sub_genomic_coord,
            sub_genomic_coord_end,
            read_samples=self.read_samples,
        )
        if compact:
            ret._compact()
        return ret
    
    def get_submatrix_from_read_mask(
        self, allowed_reads: np.ndarray(), compact: bool = True
    ) -> SparseMethylationMatrixContainer:
        """Creates a submatrix containing only the given reads and their
        methylation calls.

        :param allowed_reads: boolean array of shape (n_reads)
        :param compact: whether or not to genomic sites that no longer have
        any data
        :return: new SparseMethylationMatrixContainer object
        """
        sub_met_matrix = self.met_matrix[allowed_reads, :]
        sub_read_names = self.read_names[allowed_reads]
        sub_read_samles = self.read_samples[allowed_reads] if self.read_samples is not None else None
        
        ret = SparseMethylationMatrixContainer(
            sub_met_matrix,
            sub_read_names,
            self.genomic_coord,
            self.genomic_coord_end,
            read_samples=sub_read_samles,
        )
        if compact:
            ret._compact()
        return ret
    
    def get_submatrix_from_read_names(
        self, allowed_reads: List[str], compact: bool = True
    ) -> SparseMethylationMatrixContainer:
        """Creates a submatrix containing only the given reads and their
        methylation calls.

        :param allowed_reads: The read names to keep
        :param compact: whether or not to genomic sites that no longer have
        any data
        :return: new SparseMethylationMatrixContainer object
        """
        idx = np.array([read in allowed_reads for read in self.read_names])
        return self.get_submatrix_from_read_mask(idx, compact=compact)
    
    def get_genomic_region(self) -> Tuple[int, int]:
        """
        :return: Tuple containing min and max genomic coordinate managed by the
        matrix container object
        """
        return self.genomic_coord[0], self.genomic_coord[-1]
    
    def merge(
        self: SparseMethylationMatrixContainer,
        other: SparseMethylationMatrixContainer,
        sample_names_mode: str = "keep",
        sample_names: Union[str, Tuple[str, str]] = None,
    ):
        """
        Merge matrix with another matrix along the same genomic coordinates
        :param other: other matrix
        :param sample_names: Tuple of 2 with new sample names (see "sample_names_mode")
        :param sample_names_mode: If "keep", the sample names from both matrices are kept.
                                  If "replace" then the sample names will be replaced by the "sample_names" parameter
                                  If "append" then the sample names will be the "sample_names" parameter plus
                                  the original sample names
        :return: a new SparseMethylationMatrix
        """
        coords = sorted(list(set(self.genomic_coord).union(set(other.genomic_coord))))
        coords_self_end = [np.max(self.genomic_coord_end[self.genomic_coord == g], initial=0) for g in coords]
        coords_other = [np.max(other.genomic_coord_end[other.genomic_coord == g], initial=0) for g in coords]
        coords_end = [max(am, bm) for am, bm in zip(coords_self_end, coords_other)]
        
        coord_to_index_dict = {coords[i]: i for i in range(len(coords))}
        
        self_coo = self.met_matrix.tocoo()
        other_coo = other.met_matrix.tocoo()
        data_combined = np.concatenate((self_coo.data, other_coo.data))
        row_offset = self_coo.shape[0]
        other_row = [br + row_offset for br in other_coo.row]
        
        self_col = [coord_to_index_dict[self.genomic_coord[c]] for c in self_coo.col]
        other_col = [coord_to_index_dict[other.genomic_coord[c]] for c in other_coo.col]
        
        row_combined = np.concatenate((self_coo.row, other_row))
        col_combined = np.concatenate((self_col, other_col))
        
        met_matrix = sp.csc_matrix((data_combined, (row_combined, col_combined)))
        
        read_names = np.concatenate((self.read_names, other.read_names))
        
        if sample_names_mode == "keep":
            if self.read_samples is not None and other.read_samples is not None:
                read_samples = np.concatenate((self.read_samples, other.read_samples))
            else:
                read_samples = None
        elif sample_names_mode == "replace":
            read_samples = np.array(
                [sample_names[0]] * len(self.read_names) + [sample_names[1]] * len(other.read_names)
            )
        elif sample_names_mode == "append":
            read_samples = np.array(
                [f"{sample_names[0]}_{sn}" for sn in self.read_samples]
                + [f"{sample_names[1]}_{sn}" for sn in other.read_samples]
            )
        
        return SparseMethylationMatrixContainer(
            met_matrix, read_names, np.array(coords), np.array(coords_end), read_samples
        )
    
    def join(
        self: SparseMethylationMatrixContainer,
        other: SparseMethylationMatrixContainer,
        set_sample_based_on_sharedness=False,
    ):
        """
        HIGHLY EXPERIMENTAL: Join another matrix along the same read names. This can allow visualization
        of chimeric aliognments. The new matrix returned will be on a new - ficticious - coordinate system based
        on the genomic coordinates of the first matrix.

        The coordinates of the second matrix will be offset as such:
            y' = max(x) + (y - min(y))
        where "x" are the coordinates of the first matrix and "y" the coordinates of the second matrix

        :param other: other matrix
        :return: a new SparseMethylationMatrix
        """
        reads_only_other = set(other.read_names).difference(set(self.read_names))
        
        # new_coo = sp.coo_matrix((len(reads), len(self.genomic_coord) + len(other.genomic_coord)))
        self_coo = self.met_matrix.tocoo()
        other_coo = other.met_matrix.tocoo()
        self_read_names = list(self.read_names)
        other_read_names = list(other.read_names)
        read_names = [r for r in self_read_names] + [r for r in other_read_names if r in reads_only_other]
        read_names_other_map = {other_read_names.index(r): read_names.index(r) for r in other_read_names}
        new_other_coo_row = [read_names_other_map[r] for r in other_coo.row]
        
        new_other_coo_col = other_coo.col + np.max(self_coo.col) + 1
        
        new_other_genomic_coord = other.genomic_coord - other.genomic_coord[0] + self.genomic_coord_end[-1]
        new_other_genomic_coord_end = other.genomic_coord_end - other.genomic_coord + new_other_genomic_coord
        
        new_coo_row = np.concatenate((self_coo.row, new_other_coo_row))
        new_coo_col = np.concatenate((self_coo.col, new_other_coo_col))
        new_coo_data = np.concatenate((self_coo.data, other_coo.data))
        
        new_genomic_coord = np.concatenate((self.genomic_coord, new_other_genomic_coord))
        
        new_genomic_coord_end = np.concatenate((self.genomic_coord_end, new_other_genomic_coord_end))
        
        new_met_matrix = sp.coo_matrix((new_coo_data, (new_coo_row, new_coo_col)))
        
        sharedness = [
            "L" if r not in other_read_names else "R" if r not in self_read_names else "B" for r in read_names
        ]
        if self.read_samples is not None:
            samples = np.array(
                [s for s in self.read_samples]
                + [s for r, s in zip(other.read_names, other.read_samples) if r in reads_only_other]
            )
            if set_sample_based_on_sharedness:
                samples = np.array([f"{s}_{x}" for s, x in zip(samples, sharedness)])
        elif set_sample_based_on_sharedness:
            samples = np.array(sharedness)
        else:
            samples = None
        
        new_met_matrix = SparseMethylationMatrixContainer(
            sp.csc_matrix(new_met_matrix), read_names, new_genomic_coord, new_genomic_coord_end, read_samples=samples
        )
        
        return new_met_matrix
