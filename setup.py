from setuptools import find_packages, setup, Extension


def check_package(package = 'rpy2'):
    import importlib
    try:
        importlib.import_module(package)
    except ImportError as e:
        raise ImportError("Requires {package} to "
                        "be installed before running setup.py (pip install {package})"\
                        .format(package = package))
    finally:
        pass

_ = [check_package(p) for p in ['rpy2','numpy','pandas']]


#try:
#    import rpy2
#except ImportError:
#    raise ImportError("Requires rpy2 to "
#            "be installed before running setup.py (pip install cython)")
#try:
#    import numpy as np
#except ImportError:
#    raise ImportError("Requires numpy to "
#            "be installed before running setup.py (pip install numpy)")
#try:
#    import pandas
#except ImportError:
#    raise ImportError("Requires pandas to "
#            "be installed before running setup.py (pip install pandas)")

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
    install_requires=[
        'pandas',
        'rpy2']
)
