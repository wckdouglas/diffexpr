FROM continuumio/miniconda3:latest as base

RUN apt-get update

RUN conda config --add channels bioconda
RUN conda config --add channels default
RUN conda config --add channels anaconda
RUN conda config --add channels r
RUN conda config --add channels conda-forge

RUN conda config --set always_yes yes --set changeps1 no

RUN conda install mamba
RUN mamba install pandas tzlocal rpy2 biopython ReportLab pytest-cov codecov  bioconductor-deseq2

COPY . /opt/diffexpr

RUN Rscript /opt/diffexpr/setup.R
RUN python /opt/diffexpr/setup.py install

ENV PYTHONPATH "${PYTHONPATH}:/opt/diffexpr"
