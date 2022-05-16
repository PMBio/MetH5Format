# MetH5Format 1.1.0

[![GitHub license](https://img.shields.io/github/license/snajder-r/meth5format.svg)](https://github.com/snajder-r/meth5format/blob/master/LICENSE)
[![DOI](https://zenodo.org/badge/303672813.svg)](https://zenodo.org/badge/latestdoi/303672813)
[![Language](https://img.shields.io/badge/Language-Python3.7+-yellow.svg)](https://www.python.org/)
[![Build Status](https://travis-ci.com/snajder-r/meth5format.svg?branch=main)](https://travis-ci.com/snajder-r/meth5format)
[![Code style: black](https://img.shields.io/badge/code%20style-black-black.svg?style=flat)](https://github.com/snajder-r/black "Black (modified)")


[![PyPI version](https://badge.fury.io/py/meth5.svg)](https://badge.fury.io/py/meth5)
[![PyPI downloads](https://pepy.tech/badge/meth5)](https://pepy.tech/project/meth5)
[![Anaconda Version](https://img.shields.io/conda/v/snajder-r/meth5?color=blue)](https://anaconda.org/snajder-r/meth5)
[![Anaconda Downloads](https://anaconda.org/snajder-r/meth5/badges/downloads.svg)](https://anaconda.org/snajder-r/meth5)

MetH5 is an HDF5-based container format for methylation calls from long reads.

In the current version, the MetH5 format can store the following information:
* Log-likelihood ratio of each methylation call
* Genomic coordinates (start and end) of each methylation call
* The read name associated with each call
* Read grouping (i.e. annotation such as samples or haplotypes)

## Installation

Through pip:

```
pip install meth5
````

Through anaconda:

```
conda install -c snajder-r meth5
```

##  Usage

### Creating a MetH5 file from nanopolish methylation calls

You can create a MetH5 file with the following command, where `INPUT_PATH` refers to either a nanopolish tsv output file (may or may not be gzipped) or it can be a directory which contains only said files. 

```
meth5 create_m5 --input_paths INPUT_PATH1 [INPUT_PATH2 ...] --output_file OUTPUT_FILE.m5
```

In order to annotate reads with read grouping (for example as samples or haplotypes) you can do so by running: 

```
meth annotate_reads --m5file M5FILE.m5 --read_groups_key READ_GROUPS_KEY --read_group_file READ_GROUP_FILE
```

Where the `READ_GROUPS_KEY` is the key under which you want to store the annotation (you can store multiple read annotations), 
and `READ_GROUP_FILE` is a tab-delimited file containg read name and read group. For example:

```
read_name   group
7741f9ee-ad41-42a4-99b2-290c66960410    1
4f18b48e-a1d3-49ad-ace3-cfb96b78ad79    2
...
```

### Quick start for python API

Here an example on how to access methylation values from a MetH5 file:

```python
from meth5.meth5 import MetH5File

with MetH5File(filename, mode="r") as m:
    # List chromosomes in the MetH5 file
    m.get_chromosomes()
    
    # Access chromosome 7
    chr7 = m["chr7"]
    
    # Get number of chunks
    chr7.get_number_of_chunks()
    
    # Get a container that manages the values of chunk 3
    # (note that the data is not yet loaded into memory)
    values = chr7.get_chunk(3)
    
    # Get the log-likelihood ratios in the container as a numpy array of shape (n,)
    llrs = values.get_llrs()
    
    # Get the genomic start and end locations for each methylation call in the 
    # chunk as a numpy array  of shape (n,2) 
    ranges = values.get_ranges()
    
    # Compute methylation rate (beta-score of methylation) for each genomic location,
    # as well as the respective coordinates
    met_rates, met_rate_ranges = values.get_llr_site_rate()
    
    # You can also compute other aggregates if you like
    met_count, met_count_ranges = values.get_llr_site_aggregate(aggregation_fun=lambda llrs: (llrs>2).sum())
    
    # Instead of accessing chunk wise, you can query a genomic range
    values = chr7.get_values_in_range(36852906, 37449223)
```

A more detailed API documentation is in the works. Stay tuned!

### Sparse methylation matrix

In addition to accessing methylation calls in its unraveled form, the `meth5` library also contains a way to represent
the methylation calls as a sparse matrix. Seeing how the values are already stored in the MetH5 file in the same way a
coordinate sparse matrix would be stored in memory, this is a very cheap operation. Example:

```python
from meth5.meth5 import MetH5File

with MetH5File(filename, mode="r") as m:
    values = m["chr7"].get_values_in_range(36852906, 37449223)
    
    # The parameter "read_read_names" allows is to choose whether we want to load the actual
    # read names into memory. It's slightly more expensive than not reading it, so only load them
    # if you are interested in them
    matrix = values.to_sparse_methylation_matrix(read_read_names=True)

    # This is a scipy.sparse.csc_matrix matrix of dimension (r,s), containing the log-likelihood ratios of methylation
    # where r is the number of reads covering the genomic range we selected, and s is the number of unique genomic 
    # ranges for which we have methylation calls. Since an LLR of 0 means total uncertainty, a 0 indicates no call.
    matrix.met_matrix
    
    # A numpy array of shape (s, ) containing the start position for each unique genomic range
    matrix.genomic_coord
    # A numpy array of shape (s, ) containing the end position for each unique genomic range
    matrix.genomic_coord_end
    
    # A numpy array of shape (r, ) containing the read names
    matrix.read_names
    
    # Get a submatrix containing only the first 10 genomic locations
    submatrix = matrix.get_submatrix(0, 10)

    # Get a submatrix containing only the reads in the provided list of read names
    submatrix = matrix.get_submatrix_from_read_names(allowed_read_names)
```



## The MetH5 Format

A MetH5 file is an HDF5 container that stores methylation calls for long reads. The structure of the HDF5 file is as follows:

```
/
├─ chromosomes
│  ├─ CHROMOSOME_NAME1
│  │  ├─ llr (float dataset of shape (n,))
│  │  ├─ read_id (int dataset of shape (n,))
│  │  ├─ range (int dataset of shape (n,2))
│  │  └─ chunk_ranges (dataset of shape (c, 2))
│  │   
│  ├─ CHROMOSOME_NAME2
│  │  └─ ...
│  └─ ...
└─ reads
   ├─ read_name_mapping (string dataset of shape (r,))
   └─ read_groups
      ├─ READ_GROUP_KEY1 (int dataset of shape (r,))
      ├─ READ_GROUP_KEY2 (int dataset of shape (r,))
      └─ ... 
```

Where `n` is the number of methylation calls in the respective chromosome, `c` is the number of chunks, and `r`is the total number of reads across all chromosomes.

---

## Citing

The repository is archived at Zenodo. If you use `meth5` please cite as follow:

Rene Snajder. (2021, May 18). snajder-r/meth5. Zenodo. https://doi.org/10.5281/zenodo.4772327

## Authors and contributors

* Rene Snajder (@snajder-r): rene.snajder(at)dkfz-heidelberg.de