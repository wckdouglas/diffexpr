conda config --add channels bioconda

conda install bioconductor-deseq2 \
	bioconductor-dexseq \
	pandas \
	rpy2  \
	biopython

pip install .
