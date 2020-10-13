import setuptools

setuptools.setup(name='meth5',
      version='0.1',
      description='HDF5-based container format for Methylation calls from long reads',
      url='https://github.com/snajder-r/meth5format',
      author='Rene Snajder',
      license='MIT',
      packages=setuptools.find_packages(),
      install_requires=[
          'numpy',
          'scipy==1.4.1',
          'pandas',
          'h5py'
      ],
      python_requires='>=3.7',
      zip_safe=False)
