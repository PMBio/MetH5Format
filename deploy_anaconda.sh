#!/bin/bash

echo "compile package from setup.py"
python setup.py sdist

echo "Build noarch package for conda..."
conda-build meta.yaml --python 3.7 --output-folder conda_build -c bioconda -c conda-forge -c snajder-r --no-include-recipe

echo "Logging in to conda..."
anaconda login

echo "Deploying to conda..."
anaconda upload conda_build/**/*.tar.bz2

echo "Cleaning up"

rm -Rf dist
rm -Rf conda_build
rm -Rf build
rm -Rf *.egg-info
rm -Rf .eggs

exit 0
