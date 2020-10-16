from __future__ import annotations
from typing import List, Tuple
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

    def get_submatrix(self, start: int, end: int) -> SparseMethylationMatrixContainer:
        """Creates a submatrix containing only the genomic positions
        between start and end index.

        Note: start and end are matrix indices, not genomic coordinates
        Use get_submatrix_from_genomic_locations if you want to instead subset
        the matrix based on genomic coordinates.

        :param start: start index in matrix
        :param end: end index in matrix
        :return: new SparseMethylationMatrixContainer object
        """
        sub_met_matrix = self.met_matrix[:, start:end]
        sub_genomic_coord = self.genomic_coord[start:end]
        sub_genomic_coord_end = self.genomic_coord_end[start:end]

        ret = SparseMethylationMatrixContainer(
            sub_met_matrix, self.read_names, sub_genomic_coord, sub_genomic_coord_end, read_samples=self.read_samples,
        )
        ret._compact()
        return ret

    def get_submatrix_from_genomic_locations(self, start_base: int, end_base: int) -> SparseMethylationMatrixContainer:
        """Returns a new SparseMethylationMatrixContainer containing
        only data for the provided genomic reigon.

        :param start_base: genomic start location
        :param end_base: genomic end location
        :return: new SparseMethylationMatrixContainer object
        """
        start = self.coord_to_index_dict[start_base]
        end = self.coord_to_index_dict[end_base]
        return self.get_submatrix(start, end)

    def get_submatrix_from_read_mask(self, allowed_reads: np.ndarray()) -> SparseMethylationMatrixContainer:
        """Creates a submatrix containing only the given reads and their
        methylation calls.

        :param allowed_reads: boolean array of shape (n_reads)
        :return: new SparseMethylationMatrixContainer object
        """
        idx = allowed_reads
        sub_met_matrix = self.met_matrix[idx, :]
        sub_read_names = self.read_names[idx]
        sub_read_samles = self.read_samples[idx] if self.read_samples is not None else None

        ret = SparseMethylationMatrixContainer(
            sub_met_matrix, sub_read_names, self.genomic_coord, self.genomic_coord_end, read_samples=sub_read_samles,
        )
        ret._compact()
        return ret

    def get_submatrix_from_read_names(self, allowed_reads: List[str]) -> SparseMethylationMatrixContainer:
        """Creates a submatrix containing only the given reads and their
        methylation calls.

        :param allowed_reads: The read names to keep
        :return: new SparseMethylationMatrixContainer object
        """
        idx = [read in allowed_reads for read in self.read_names]
        return self.get_submatrix_from_read_mask(idx)

    def get_genomic_region(self) -> Tuple[int, int]:
        """
        :return: Tuple containing min and max genomic coordinate managed by the
        matrix container object
        """
        return self.genomic_coord[0], self.genomic_coord[-1]
