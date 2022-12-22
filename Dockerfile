FROM python:3.9.16-bullseye AS base


RUN apt-get update \
    && apt-get install -y r-base r-base-dev r-bioc-deseq2 libcurl4-openssl-dev libxml2-dev libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


RUN pip install pandas tzlocal \
    biopython ReportLab pytest rpy2

# R packages
COPY setup.R /opt/setup.R
RUN /usr/bin/Rscript /opt/setup.R

FROM base AS diffexpr

COPY . /opt/diffexpr

RUN python /opt/diffexpr/setup.py install

ENV PYTHONPATH "${PYTHONPATH}:/opt/diffexpr"
WORKDIR /opt/diffexpr
RUN pytest -vvv 

FROM diffexpr AS diffexpr_dev
RUN pip install jupyterlab matplotlib seaborn
RUN /usr/bin/Rscript -e "install.packages('IRkernel')"
WORKDIR /opt/diffexpr
RUN pytest -vvv 

CMD ["jupyter-lab", "--allow-root", "--ip", "0.0.0.0", "--port", "1234"]
