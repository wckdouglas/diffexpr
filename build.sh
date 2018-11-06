conda config --add channels r
conda config --add channels bioconda
conda config --add channels auto
conda config --add channels conda-forge

conda install bioconductor-deseq2 \
	bioconductor-dexseq \
	pandas \
	rpy2  \
	biopython

pip install .