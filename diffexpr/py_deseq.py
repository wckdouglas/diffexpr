import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, numpy2ri, Formula
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr
import numpy as np
import logging
logging.basicConfig(level = logging.INFO)
logger = logging.getLogger('DESeq2')
deseq = importr('DESeq2')
'''
Adopted from: https://stackoverflow.com/questions/41821100/running-deseq2-through-rpy2
'''

to_dataframe = robjects.r('function(x) data.frame(x)')

class py_DESeq2:
    '''
    DESeq2 object through rpy2

    Args:
        count_matrix (pd.DataFrame): should be a pandas dataframe with each column as count, and a id column for gene id
        design_matrix (pd.DataFrame): an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
        design_formula (str): see DESeq2 manual, example: "~ treatment""
        gene_column (str): column name of gene id columns, example "id"


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

    '''
    def __init__(self, count_matrix, design_matrix, design_formula, gene_column='id'):
        try:
            assert gene_column in count_matrix.columns, 'Wrong gene id column name'
            gene_id = count_matrix[gene_column]
        except AttributeError:
            sys.exit('Wrong Pandas dataframe?')

        self.dds = None
        self.result = None
        self.deseq_result = None
        self.resLFC = None
        self.comparison = None
        self.normalized_count_df = None
        self.gene_column = gene_column
        self.gene_id = count_matrix[self.gene_column]
        self.samplenames = count_matrix.columns[count_matrix.columns != self.gene_column]
        with localconverter(robjects.default_converter + pandas2ri.converter):
            self.count_matrix = robjects.conversion.py2rpy(count_matrix.set_index(self.gene_column))
            self.design_matrix = robjects.conversion.py2rpy(design_matrix)
        self.design_formula = Formula(design_formula)
        self.dds = deseq.DESeqDataSetFromMatrix(countData=self.count_matrix, 
                                        colData=self.design_matrix,
                                        design=self.design_formula)


    def run_deseq(self, **kwargs):
        """
        actually running deseq2

        Args:
            **kwargs: Any keyword arguments for DESeq

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
        """
        
        for key, value in kwargs.items():
            if key == 'reduced':
                kwargs[key] = Formula(value)
        self.dds = deseq.DESeq(self.dds, **kwargs)
        self.comparison = list(deseq.resultsNames(self.dds))


    def get_deseq_result(self, contrast=None, **kwargs):
        '''
        DESeq2: result(dds, contrast)
        making a dds.deseq_result pandas dataframe
        '''

        if contrast:
            if len(contrast)==3:
                R_contrast = robjects.vectors.StrVector(np.array(contrast)) 
            else:
                if len(contrast) != 2:
                    raise ValueError('Contrast must be length of 3 or 2')
                R_contrast = robjects.ListVector({None:con for con in contrast})
            logger.info('Using contrast: %s' %contrast)
            self.result = deseq.results(self.dds, contrast = R_contrast, **kwargs) # Robject
        else:
            self.result = deseq.results(self.dds, **kwargs) # R object
        self.deseq_result = to_dataframe(self.result) # R dataframe
        with localconverter(robjects.default_converter + pandas2ri.converter):
            self.deseq_result = robjects.conversion.rpy2py(self.deseq_result) ## back to pandas dataframe
        self.deseq_result[self.gene_column] = self.gene_id.values

    def normalized_count(self):
        '''
        Returns:
            pd.DataFrame: a dataframe in the format of DESeq2::Counts(dds, normalized=TRUE) in R
        '''
        normalized_count_matrix = deseq.counts_DESeqDataSet(self.dds, normalized=True)
        normalized_count_matrix = to_dataframe(normalized_count_matrix)
        # switch back to python
        with localconverter(robjects.default_converter + pandas2ri.converter):
            self.normalized_count_df = robjects.conversion.rpy2py(normalized_count_matrix)  
        self.normalized_count_df[self.gene_column] = self.gene_id.values
        logger.info('Normalizing counts')
        return self.normalized_count_df

    def lfcShrink(self, coef, method = 'apeglm'):
        '''
        Perform LFC shrinkage on the DDS object
        see: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

        Be sure to check dds.comparison to see which coef (1-base) to use

        Args:
            coef (int): 1-based index for selecting which of dds.comparison to show
            method (str): DESeq2 lfcshrink method ("apeglm", "ashr", "normal")

        Returns:
            pandas.DataFrame: a deseq2 result table
        '''
        lfc = deseq.lfcShrink(self.dds, res = self.result, coef = coef, type = method)
        with localconverter(robjects.default_converter + pandas2ri.converter):
            lfc = robjects.conversion.rpy2py(to_dataframe( lfc))

        return lfc\
            .reset_index()\
            .rename(columns = {'index':self.gene_column})
    
