FROM continuumio/miniconda3:latest AS base

RUN apt-get update

RUN conda config --add channels bioconda
RUN conda config --add channels default
RUN conda config --add channels anaconda
RUN conda config --add channels r
RUN conda config --add channels conda-forge

RUN conda config --set always_yes yes --set changeps1 no

RUN conda install mamba
RUN mamba install pandas tzlocal \
  rpy2 biopython ReportLab \
  bioconductor-deseq2 

FROM base AS diffexpr

COPY . /opt/diffexpr

RUN Rscript /opt/diffexpr/setup.R
RUN python /opt/diffexpr/setup.py install
RUN conda clean --all --yes

ENV PYTHONPATH "${PYTHONPATH}:/opt/diffexpr"

FROM diffexpr AS diffexpr_dev
RUN mamba install jupyterlab matplotlib seaborn
RUN conda clean --all --yes
CMD ["jupyter-lab", "--allow-root", "--ip", "0.0.0.0", "--port", "1234"]
