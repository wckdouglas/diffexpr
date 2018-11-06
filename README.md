# diffexpr # 
[![Build Status](https://travis-ci.org/wckdouglas/diffexpr.svg?branch=master)](https://travis-ci.org/wckdouglas/diffexpr)

A python package using ```rpy2``` to port [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and [DEXSeq](https://bioconductor.org/packages/release/bioc/html/DEXSeq.html) into python.

## INSTALL ##
Dependencies are ```pandas``` (python), ```rpy2``` (python), ```DESeq2``` (R) and ```DEXSeq``` (R)
Best way to build dependencies should be via conda. 

```
conda config --add channels r
conda config --add channels bioconda
conda config --add channels auto
conda config --add channels conda-forge

conda install bioconductor-deseq2 \
    bioconductor-dexseq \
    pandas \
    rpy2 
```

## Example ##
An example of running DESeq2 in *python* using ```diffexp``` package is provided [here](https://github.com/wckdouglas/diffexp/blob/master/example/deseq_example.ipynb).
