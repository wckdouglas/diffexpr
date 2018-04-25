
from __future__ import print_function
import pandas as pd
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, Formula
pandas2ri.activate()
from rpy2.robjects.packages import importr
dexseq = importr('DEXSeq')
bp = importr('BiocParallel')
'''
Adopted from: https://stackoverflow.com/questions/41821100/running-deseq2-through-rpy2
'''

to_dataframe = robjects.r('function(x) data.frame(x)')

class py_DEXSeq:
    '''
    DEXSeq2 object through rpy2

    input:
    count_matrix: should be a pandas dataframe with each column as count, and a id column for exon id
        example:
        exonID    sampleA    sampleB
        geneA    5    1
        geneB    4    5
        geneC    1    2
    design_matrix: an design matrix in the form of pandas dataframe, see DESeq2 manual, samplenames as rownames
                treatment
    sampleA1        A
    sampleA2        A
    sampleB1        B
    sampleB2        B

    design_formula: see DEXSeq manual, example: "~ sample + exon + exon:treatment""
    var_column: will pass to fitExpToVar for DEXSeq exon fold change
    exons: exon id
    genes: gene id for dexseq grouping
    threads: number of threads to use
    '''
    def __init__(self, count_matrix, design_matrix, design_formula, 
                var_column = 'condition',
                exons=None, genes=None, threads=1):
        try:
            assert var_column in design_matrix.columns, 'Wrong var column for DEXSeq'
        except AttributeError:
            sys.exit('Wrong Pandas dataframe?')

        self.dxd = None
        self.dxd_res = None
        self.dexseq_result = None
        self.comparison = None
        self.normalized_count_matrix = None
        self.exons = exons
        self.genes = genes
        self.sample_names = count_matrix.columns
        self.count_matrix = pandas2ri.py2ri(count_matrix)
        self.design_matrix = pandas2ri.py2ri(design_matrix)
        self.design_formula = Formula(design_formula)
        self.BPPARAM = bp.MulticoreParam(workers=threads)
        self.var_column = var_column


    def run_dexseq(self, **kwargs):

        self.dxd = dexseq.DEXSeqDataSet(countData=self.count_matrix, 
                                        sampleData=self.design_matrix,
                                        design=self.design_formula,
                                        featureID = self.exons,
                                        groupID = self.genes)
        print('Constructed DXD object')
        self.dxd = dexseq.estimateSizeFactors_DEXSeqDataSet(self.dxd) 
        self.dxd = dexseq.estimateDispersions_DEXSeqDataSet(self.dxd, BPPARAM=self.BPPARAM)
        print('Starting DEXSeq test')
        self.dxd = dexseq.testForDEU(self.dxd, BPPARAM=self.BPPARAM)
        self.dxd = dexseq.estimateExonFoldChanges(self.dxd, 
                                                fitExpToVar=self.var_column,
                                                BPPARAM=self.BPPARAM) 
        print('Finished DEXSeq fold change')

    def get_dexseq_result(self, **kwargs):
        self.deu_res = dexseq.DEXSeqResults(self.dxd)
        self.dexseq_result = to_dataframe(self.deu_res, **kwargs)
        self.dexseq_result = pandas2ri.ri2py(self.dexseq_result) ## back to pandas dataframe 
        self.dexseq_result = self.dexseq_result.reset_index()
        self.dexseq_result.drop('genomicData', axis=1, inplace=True)

    def normalized_count(self):
        self.normalized_count_matrix = dexseq.counts_DEXSeqResults(self.deu_res,normalized=True)
        self.normalized_count_matrix = pandas2ri.ri2py(self.normalized_count_matrix)
        self.normalized_count_matrix = pd.DataFrame(self.normalized_count_matrix,
                                                    columns = self.sample_names) \
                .rename(columns = lambda x: 'normCount_'+x)
        return self.normalized_count_matrix
