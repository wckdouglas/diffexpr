"""
Running DESeq2 from python via rpy2

Adopted from: https://stackoverflow.com/questions/41821100/running-deseq2-through-rpy2

If there are any functions that is missing for your need, feel free to file an 
[issue](https://github.com/wckdouglas/diffexpr/issues) or even better, make a 
[PR](https://github.com/wckdouglas/diffexpr/pulls)!

A jupyter notebook with examples can be found at:
https://github.com/wckdouglas/diffexpr/blob/master/example/deseq_example.ipynb
"""

import logging
from typing import Dict

import numpy as np
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import Formula, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

# setup logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("DESeq2")

# R packages as python objects
r_utils = importr("utils")
deseq = importr("DESeq2")
tximport = importr("tximport")
multicore = importr("BiocParallel")
summarized_experiment = importr("SummarizedExperiment")

# get version of deseq2
_DESEQ2_VERSION_INT = r_utils.packageVersion("DESeq2")
DESEQ2_VERSION = ".".join(map(str, robjects.conversion.rpy2py(_DESEQ2_VERSION_INT)[0]))

# a R function to make matrix into dataframe
to_dataframe = robjects.r("function(x) data.frame(x)")


class py_DESeq2:
    """
    DESeq2 object through rpy2

    Args:
        count_matrix (Union[pd.DataFrame, Dict[str,str]): should be a pandas dataframe with each column as count, and a id column for gene id,
            unless kallisto=True, then this is expected to be a dictionary of key: sample name, value: abundance.h5 file
        design_matrix (pd.DataFrame): an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
        design_formula (str): see DESeq2 manual, example: "~ treatment""
        gene_column (str): column name of gene id columns (default: "id")
        threads (int): how many threads to used in running deseq, if threads > 1 is provided,
            `parallel=True` will be used in `DESeq2::DESeq`, `DESeq2::results`, and `DESeq2::lfcShrink`
            (default: 1)


    count_matrix example::

        id    sampleA    sampleB
        geneA    5    1
        geneB    4    5
        geneC    1    2

    Design matrix example::

                    treatment
        sampleA1        A
        sampleA2        A
        sampleB1        B
        sampleB2        B

    """

    def __init__(
        self, count_matrix, design_matrix, design_formula, gene_column="id", threads=1, kallisto=False, tx2gene=None
    ):
        if not isinstance(threads, int):
            raise ValueError("threads must be an integer")
        multicore.register(multicore.MulticoreParam(threads))

        # set up the deseq2 object
        self.dds = None
        self.result = None
        self.deseq_result = None
        self.resLFC = None
        self.comparison = None
        self.normalized_count_df = None
        self.parallel = threads > 1
        self.gene_id = None
        self.gene_column = None

        if kallisto:
            if tx2gene is None:
                raise ValueError("tx2gene must be specified")
            self.from_kallisto(count_matrix, design_matrix, design_formula, tx2gene)
        else:
            self.init_matrix(count_matrix, design_matrix, design_formula, gene_column)

    def init_matrix(self, count_matrix, design_matrix, design_formula, gene_column):
        """
        Initialize deseq from count matrix

        Args:
            count_matrix (pd.DataFrame):
            design_matrix (pd.DataFrame):
            design_formula (str):
            gene_column (str):
        """
        # input validation
        for df in [count_matrix, design_matrix]:
            if not isinstance(df, pd.DataFrame):
                raise ValueError("count_matrix and design_matrix should be pd.DataFrame type")

        if gene_column not in count_matrix.columns:
            raise ValueError("The given gene_column name is not a column in  count_matrix dataframe")

        self.gene_column = gene_column
        self.gene_id = count_matrix[self.gene_column]
        self.samplenames = count_matrix.columns[count_matrix.columns != self.gene_column]
        with localconverter(robjects.default_converter + pandas2ri.converter):
            self.count_matrix = robjects.conversion.py2rpy(count_matrix.set_index(self.gene_column))
            self.design_matrix = robjects.conversion.py2rpy(design_matrix)
        self.design_formula = Formula(design_formula)
        self.dds = deseq.DESeqDataSetFromMatrix(
            countData=self.count_matrix, colData=self.design_matrix, design=self.design_formula
        )

    def from_kallisto(
        self, h5_file_list: Dict[str, str], design_matrix: pd.DataFrame, design_formula: str, tx2gene: pd.DataFrame
    ):
        """
        Initialize deseq from Tximport kallisto files

        :param h5_file_list: dictionary of key: sample name, value: abundance.h5 file
        :param design_matrix: an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
        :param str design_formula: see DESeq2 manual, example: "~ treatment""
        """
        files = robjects.StrVector(list(h5_file_list.values()))
        files.names = list(h5_file_list.keys())
        self.design_formula = Formula(design_formula)
        with localconverter(robjects.default_converter + pandas2ri.converter):
            self.design_matrix = robjects.conversion.py2rpy(design_matrix)
            tx2gene = robjects.conversion.py2rpy(tx2gene)
        self.txi = tximport.tximport(
            files, type="kallisto", txOut=False, tx2gene=tx2gene, countsFromAbundance="scaledTPM"
        )
        logger.info(f"Read kallisto files: {files}")
        self.dds = deseq.DESeqDataSetFromTximport(self.txi, colData=self.design_matrix, design=self.design_formula)

    def run_deseq(self, **kwargs):
        """
        actually running deseq2 and setup the dds and comparison
        fields in the object

        From DESeq2 manual:

        DESeq(
            object,
            test = c("Wald", "LRT"),
            fitType = c("parametric", "local", "mean", "glmGamPoi"),
            sfType = c("ratio", "poscounts", "iterate"),
            betaPrior,
            full = design(object),
            reduced,
            quiet = FALSE,
            minReplicatesForReplace = 7,
            modelMatrixType,
            useT = FALSE,
            minmu = if (fitType == "glmGamPoi") 1e-06 else 0.5,
            parallel = FALSE,
            BPPARAM = bpparam()
        )

        Args:
            **kwargs: Any keyword arguments for DESeq
        Returns:
            NoneType

        """

        for key, value in kwargs.items():
            if key == "reduced":
                kwargs[key] = Formula(value)
            if key == "parallel":
                raise ValueError("parallel is inferred from the provided thread count")
        self.dds = deseq.DESeq(self.dds, parallel=self.parallel, **kwargs)
        self.comparison = list(deseq.resultsNames(self.dds))

    def get_deseq_result(self, contrast=None, **kwargs):
        """
        DESeq2: result(dds, contrast)
        making a dds.deseq_result pandas dataframe

        Args:
            contrast (list[str]): list of string annotating the contrast to compute
                (see DESeq2 manual http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#contrasts)
            **kwargs (Dict[str,Any]): any other parameters to pass into DESeq2::results function
        Returns:
            NoneType
        """

        if contrast:
            if len(contrast) == 3:
                r_contrast = robjects.vectors.StrVector(np.array(contrast))
            else:
                if len(contrast) != 2:
                    raise ValueError("Contrast must be length of 3 or 2")
                r_contrast = robjects.ListVector({None: con for con in contrast})
            logger.info("Using contrast: %s" % contrast)
            self.result = deseq.results(self.dds, contrast=r_contrast, parallel=self.parallel, **kwargs)  # Robject
        else:
            self.result = deseq.results(self.dds, **kwargs)  # R object
        self.deseq_result = to_dataframe(self.result)  # R dataframe
        with localconverter(robjects.default_converter + pandas2ri.converter):
            self.deseq_result = robjects.conversion.rpy2py(self.deseq_result)  ## back to pandas dataframe

        if self.gene_column is not None:
            self.deseq_result[self.gene_column] = self.gene_id.values

    def normalized_count(self):
        """
        Returns a normalized count data frame

        Returns:
            pd.DataFrame: a dataframe in the format of DESeq2::Counts(dds, normalized=TRUE) in R
        """
        normalized_count_matrix = deseq.counts_DESeqDataSet(self.dds, normalized=True)
        normalized_count_matrix = to_dataframe(normalized_count_matrix)
        # switch back to python
        with localconverter(robjects.default_converter + pandas2ri.converter):
            self.normalized_count_df = robjects.conversion.rpy2py(normalized_count_matrix)

        if self.gene_column is not None:
            self.normalized_count_df[self.gene_column] = self.gene_id.values
        logger.info("Normalizing counts")
        return self.normalized_count_df

    def lfcShrink(self, coef, method="apeglm", **kwargs):
        """
        Perform LFC shrinkage on the DDS object
        see: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

        Be sure to check dds.comparison to see which coef (1-base because it's passing into the R code)
        to use

        Args:
            coef (int): 1-based index for selecting which of dds.comparison to show
            method (str): DESeq2 lfcshrink method ("apeglm", "ashr", "normal")

        Returns:
            pandas.DataFrame: a deseq2 result table
        """
        lfc = deseq.lfcShrink(self.dds, res=self.result, coef=coef, type=method, parallel=self.parallel, **kwargs)
        with localconverter(robjects.default_converter + pandas2ri.converter):
            lfc = robjects.conversion.rpy2py(to_dataframe(lfc))

        return lfc.reset_index().rename(columns={"index": self.gene_column})

    def vst(self, blind=True, fit_type="parametric"):
        """
        deseq varianceStabilizingTransformation
        see: https://rdrr.io/bioc/DESeq2/man/varianceStabilizingTransformation.html

        essentially running R code:

        >>> vsd = DESeq2::varianceStabilizingTransformation(dds, blind=True, fitType="parametric")
        >>> SummarizedExperiment::assay(vsd)

        Example:

        >>> dds = py_DESeq2(
                count_matrix=df,
                design_matrix=sample_df,
                design_formula="~ batch + sample",
                gene_column="id",
        )
        >>> dds.vst(blind=True, fit_type="parametric")


        Args:
            blind (bool):  whether to blind the transformation to the experimental design
            fit_type (str): should be either "parametric", "local", "mean"
        Returns:
            pandas.DataFrame: a vst transformed count table
        """
        if self.dds is None:
            raise ValueError("Empty DESeq object")

        acceptable_fit_types = set(["parametric", "local", "mean"])
        if fit_type not in acceptable_fit_types:
            raise ValueError(f"fit_type must be {acceptable_fit_types}")

        vst_matrix = summarized_experiment.assay(
            deseq.varianceStabilizingTransformation(self.dds, blind=blind, fitType=fit_type)
        )
        vst_df = to_dataframe(vst_matrix)
        with localconverter(robjects.default_converter + pandas2ri.converter):
            vst_counts = robjects.conversion.rpy2py(vst_df)

        logger.info("Processed variance stablizing transformation")
        return vst_counts.reset_index().rename(columns={"index": self.gene_column})

    def rlog(self, blind=True, fit_type="parametric"):
        """
        deseq rlog
        see: https://rdrr.io/bioc/DESeq2/man/rlog.html

        TODO: DESeq2 version of this function accepts two additional optional arguments
        'intercept' and 'betaPriorVar' that have not been explicitly ported here.

        essentially running R code:

        >>> rld = DESeq2::rlog(dds, blind=True, fitType="parametric")
        >>> SummarizedExperiment::assay(rld)

        Example:

        >>> dds = py_DESeq2(
                count_matrix=df,
                design_matrix=sample_df,
                design_formula="~ batch + sample",
                gene_column="id",
        )
        >>> dds.rlog(blind=True, fit_type="parametric")


        Args:
            blind (bool):  whether to blind the transformation to the experimental design
            fit_type (str): should be either "parametric", "local", "mean"
        Returns:
            pandas.DataFrame: a rlog transformed count table
        """
        if self.dds is None:
            raise ValueError("Empty DESeq object")

        acceptable_fit_types = set(["parametric", "local", "mean"])
        if fit_type not in acceptable_fit_types:
            raise ValueError(f"fit_type must be {acceptable_fit_types}")

        rlog_matrix = summarized_experiment.assay(deseq.rlog(self.dds, blind=blind, fitType=fit_type))
        rlog_df = to_dataframe(rlog_matrix)
        with localconverter(robjects.default_converter + pandas2ri.converter):
            rlog_counts = robjects.conversion.rpy2py(rlog_df)

        logger.info("Processed rlog transformation")
        return rlog_counts.reset_index().rename(columns={"index": self.gene_column})

    @property
    def deseq2_version(self):
        """
        Return DESeq2 version number

        :return: string
        """
        return DESEQ2_VERSION
