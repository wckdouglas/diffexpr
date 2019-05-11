#!/usr/bin/env python
import os
import pandas as pd 
import numpy as np
from diffexpr.py_deseq import py_DESeq2
import warnings
warnings.filterwarnings("ignore")
test_data_path = os.path.dirname(os.path.realpath(__file__)) + '/data'


def test_deseq():

    df = pd.read_csv(test_data_path + '/ercc.tsv', sep='\t')
    """
        id     A_1     A_2     A_3     B_1     B_2     B_3
    0  ERCC-00002  111461  106261  107547  333944  199252  186947
    1  ERCC-00003    6735    5387    5265   13937    8584    8596
    2  ERCC-00004   17673   13983   15462    5065    3222    3353
    3  ERCC-00009    4669    4431    4211    6939    4155    3647
    4  ERCC-00012       0       2       0       0       0       0
    """


    sample_df = pd.DataFrame({'samplename': df.columns}) \
        .query('samplename != "id"')\
        .assign(sample = lambda d: d.samplename.str.extract('([AB])_', expand=False)) \
        .assign(batch = lambda d: d.samplename.str.extract('_([123])', expand=False)) 
    sample_df.index = sample_df.samplename

    dds = py_DESeq2(count_matrix = df,
               design_matrix = sample_df,
               design_formula = '~ sample',
               gene_column = 'id')
    
    dds.run_deseq() 
    dds.get_deseq_result()
    res = dds.deseq_result 
    res.query('padj < 0.05')
    assert(res.query('padj < 0.05').shape == (43,7))

    res.to_csv(test_data_path + '/py_deseq.tsv', index=False, sep='\t')


def test_result():
    os.chdir(os.path.dirname(test_data_path))
    os.system('Rscript deseq.R')

    py = pd.read_table(test_data_path + '/py_deseq.tsv') 
    R = pd.read_table(test_data_path + '/R_deseq.tsv') 

    for col in py.columns:
        if py.columns.dtype == 'float64':
            assert(np.all(np.isclose(py[col].fillna(0), R[col].fillna(0))))
