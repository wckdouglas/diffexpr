FROM continuumio/miniconda3:latest AS base


RUN apt-get update \
    && apt-get install -y r-base r-base-dev r-bioc-deseq2 libcurl4-openssl-dev libxml2-dev libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


# python package (diffexpr)
RUN conda config --add channels bioconda
RUN conda config --add channels default
RUN conda config --add channels anaconda
RUN conda config --add channels conda-forge
RUN conda config --set always_yes yes --set changeps1 no

RUN conda install mamba
RUN mamba install pandas tzlocal \
    biopython ReportLab pytest
RUN pip install rpy2

# R packages
COPY setup.R /opt/setup.R
RUN /usr/bin/Rscript /opt/setup.R

FROM base AS diffexpr

COPY . /opt/diffexpr

RUN python /opt/diffexpr/setup.py install
RUN conda clean --all --yes

ENV PYTHONPATH "${PYTHONPATH}:/opt/diffexpr"
WORKDIR /opt/diffexpr
RUN pytest -vvv 

FROM diffexpr AS diffexpr_dev
RUN mamba install jupyterlab matplotlib seaborn
RUN conda clean --all --yes
RUN RScript.exe  -e "install.packages('IRkernel')"
WORKDIR /opt/diffexpr
RUN pytest -vvv 

CMD ["jupyter-lab", "--allow-root", "--ip", "0.0.0.0", "--port", "1234"]
