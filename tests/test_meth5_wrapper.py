import unittest
from pathlib import Path
from tempfile import TemporaryFile
from meth5.meth5_wrapper import MetH5File, create_sparse_matrix_from_samples


class TestMetH5Format(unittest.TestCase):
    def __init__(self, *args):
        super().__init__(*args)
        self.datadir = Path("../testdata/")
        # Just some arbitrary grouping of reads
        self.read_groups = {
            "4f18b48e-a1d3-49ad-ace3-cfb96b79ad79": 1,
            "6fc0ec64-9919-4a87-8431-76d06b151888": 2,
            "7741f9ee-ad41-42a4-99b2-290c66960410": 1,
            "efba773d-fcff-48c5-a833-e3f38cafcedf": 2,
        }

    def test_create_h5(self):
        with TemporaryFile() as tmp_f:
            #  ==== Creating a new H5 file from a regular nanopolish output ====
            # Note that chunk size is only that small for the sake of the test
            with MetH5File(tmp_f, "w", chunk_size=10) as mf:
                mf.parse_and_add_nanopolish_file(
                    self.datadir.joinpath("nanopolish_calls.tsv.gz")
                )
                # Creating this index will make random access MUCH faster for large
                # H5 files
                mf.create_chunk_index()
                # Test adding read groups
                mf.annotate_read_groups("test_group", self.read_groups)

            # ==== Test if we can read the new H5 file we just created ====
            with MetH5File(tmp_f, "r", chunk_size=10) as mf:
                chroms = mf.get_chromosomes()

                # Our test data only contains a single chromosome (chr8)
                self.assertEqual(len(chroms), 1)
                self.assertEqual(chroms[0], "8")

                # Access values for that chromosome
                chrom_container = mf[chroms[0]]

                # ---- Access chunk-wise: ----
                value_container = chrom_container.get_chunk(0)
                self.assertEqual(len(value_container.get_read_groups("test_group")), 10)

                # Test conversion to sparse matrix
                matrix_container = value_container.to_sparse_methylation_matrix(
                    read_groups_key="test_group"
                )
                # Only 1 read in that chunk
                self.assertEqual(matrix_container.met_matrix.shape[0], 1)
                # 10 sites in that chunk
                self.assertEqual(matrix_container.met_matrix.shape[1], 10)

                # ---- Access region-wise: ----
                value_container = chrom_container.get_values_in_range(
                    97732352, 97745195
                )

                # Test conversion to sparse matrix
                matrix_container = value_container.to_sparse_methylation_matrix(
                    read_groups_key="test_group"
                )

                # All 4 reads in that region
                self.assertEqual(matrix_container.met_matrix.shape[0], 4)
                # 184 sites in that region
                self.assertEqual(matrix_container.met_matrix.shape[1], 184)

                # Test subsetting by extracting reads from a sample:
                sub_matrix = matrix_container.get_submatrix_from_read_mask(
                    matrix_container.read_samples == 1
                )

                # Only 2 of the reads passed the filter
                self.assertEqual(sub_matrix.met_matrix.shape[0], 2)
                # Subsetting reads also compacted the sites dimension:
                self.assertEqual(sub_matrix.met_matrix.shape[1], 146)

                # Creates a matrix from two value container (the same one - just for
                # testing)
                two_sample_matrix = create_sparse_matrix_from_samples(
                    {"A": value_container, "B": value_container},
                    sample_prefix_readnames=True,
                )
                # Matrix should have twice as many reads in two samples but the same
                # number of sites
                self.assertEqual(len(two_sample_matrix.read_samples), 8)
                self.assertEqual(len(set(two_sample_matrix.read_samples)), 2)
                self.assertEqual(two_sample_matrix.met_matrix.shape[1], 184)


if __name__ == "__main__":
    unittest.main()
