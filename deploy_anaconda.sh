#!/bin/bash
# -*- coding: utf-8 -*-
set -e

echo "Set up conda package manager"
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh --quiet
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no --set anaconda_upload no
conda update -q conda

echo "Build noarch package for conda..."
conda-build meta.yaml --python 3.7 --output-folder conda_build -c bioconda -c conda-forge --no-include-recipe

echo "compile package from setup.py"
python setup.py sdist

echo "Deploying to conda..."
anaconda -v -t $1 upload conda_build/**/*.tar.bz2

echo "Cleaning up"

exit 0
