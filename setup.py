from setuptools import find_packages, setup, Extension
#from distutils.core import setup, Extension
import glob

try:
    import rpy2
except ImportError:
    raise ImportError("Requires rpy2 to "
            "be installed before running setup.py (pip install cython)")
try:
    import numpy as np
except ImportError:
    raise ImportError("Requires numpy to "
            "be installed before running setup.py (pip install numpy)")
try:
    import pandas
except ImportError:
    raise ImportError("Requires pysam to "
            "be installed before running setup.py (pip install pysam)")

setup(
    name = 'diffexp',
    version = '0.1',
    description = 'Tools for porting differential expression DESeq to python',
    url = '',
    author = 'Douglas C. Wu',
    author_email = 'wckdouglas@gmail.com',
    license = 'MIT',
    packages = find_packages(),
    zip_safe = False,
    ext_modules = ext_modules,
    cmdclass = {'build_ext': build_ext}
)
