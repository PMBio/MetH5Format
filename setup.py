#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Long description from README file
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect info in a dictionary for setup.py
setup(
    name="meth5",
    description="HDF5 based file format for storage, retrieval, and analysis of modification predictions from Nanopore",
    version="1.0.4",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/snajder-r/meth5format",
    author="Rene Snajder",
    author_email="r.snajder@dkfz-heidelberg.de",
    license="MIT",
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3"
    ],
    install_requires=[
        "numpy>=1.19.2",
        "scipy==1.4.1",
        "pandas>=1.1.3",
        "h5py>=2.10.0",
        "h5py<3.3.0",
        "tqdm",
        "setuptools"
    ],
    packages=["meth5"],
    package_dir={"meth5": "meth5"},
    entry_points={"console_scripts": ["meth5=meth5.__main__:main"]},
)
