FROM python:3.9.16-bullseye AS base

# installation of R and associated packages (DESeq2)
RUN apt-get update \
    && apt-get install -y r-base r-base-dev r-bioc-deseq2 libcurl4-openssl-dev libxml2-dev libssl-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# installation of python dependencies
RUN pip install pandas tzlocal \
    biopython ReportLab pytest rpy2


# At this stage, we will build the diffexpr package
# and any R packages that are not installed by apt-get 
# e.g. apeglm
FROM base AS diffexpr

COPY . /opt/diffexpr
RUN python /opt/diffexpr/setup.py install
RUN /usr/bin/Rscript -e "BiocManager::install(c('apeglm'))"

# make sure python kernel knows about the diffexpr package
# TODO: do we really need this?
ENV PYTHONPATH "${PYTHONPATH}:/opt/diffexpr"

# run a test to make sure things are installed correctly
WORKDIR /opt/diffexpr
RUN pytest -vvv  

# for the dev stage, we will install packages that help
# differential expression analyses
FROM diffexpr AS diffexpr_dev

# jupyter lab and plotting
RUN pip install jupyterlab matplotlib seaborn

# R kernel for debugging
RUN /usr/bin/Rscript -e "install.packages('IRkernel'); IRkernel::installspec(user = FALSE)"

# again, run a test to make sure things at this stage are installed correctly
WORKDIR /opt/diffexpr
RUN pytest -vvv 

# docker run will spin up a jupyter lab instance at port 1234
# docker run -p 1234:1234 should be allow access from outside the 
# container
CMD ["jupyter-lab", "--allow-root", "--ip", "0.0.0.0", "--port", "1234"]
