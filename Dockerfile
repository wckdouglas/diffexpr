FROM continuumio/miniconda3:latest AS base

RUN apt-get update \
    && apt-get install -y r-base libcurl4-openssl-dev libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY setup.R /opt/setup.R

RUN Rscript /opt/setup.R

RUN conda config --add channels bioconda
RUN conda config --add channels default
RUN conda config --add channels anaconda
RUN conda config --add channels conda-forge

RUN conda config --set always_yes yes --set changeps1 no

RUN conda install mamba
RUN mamba install pandas tzlocal \
  rpy2 biopython ReportLab 

FROM base AS diffexpr

COPY . /opt/diffexpr

RUN python /opt/diffexpr/setup.py install
RUN conda clean --all --yes

ENV PYTHONPATH "${PYTHONPATH}:/opt/diffexpr"
WORKDIR /opt/diffexpr
RUN pytest -vvv

FROM diffexpr AS diffexpr_dev
RUN mamba install jupyterlab matplotlib seaborn r-recommended r-irkernel
RUN conda clean --all --yes
CMD ["jupyter-lab", "--allow-root", "--ip", "0.0.0.0", "--port", "1234"]
