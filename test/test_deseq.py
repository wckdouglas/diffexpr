#!/usr/bin/env python
import os
import warnings

import numpy as np
import pandas as pd
import pytest
from diffexpr.py_deseq import py_DESeq2

warnings.filterwarnings("ignore")
test_data_path = os.path.dirname(os.path.realpath(__file__)) + "/data"


@pytest.fixture(scope="module")
def run_r():
    os.chdir(os.path.dirname(test_data_path))
    os.system("Rscript deseq.R")


def test_deseq():

    df = pd.read_csv(test_data_path + "/ercc.tsv", sep="\t")
    """
        id     A_1     A_2     A_3     B_1     B_2     B_3
    0  ERCC-00002  111461  106261  107547  333944  199252  186947
    1  ERCC-00003    6735    5387    5265   13937    8584    8596
    2  ERCC-00004   17673   13983   15462    5065    3222    3353
    3  ERCC-00009    4669    4431    4211    6939    4155    3647
    4  ERCC-00012       0       2       0       0       0       0
    """

    sample_df = (
        pd.DataFrame({"samplename": df.columns})
        .query('samplename != "id"')
        .assign(sample=lambda d: d.samplename.str.extract("([AB])_", expand=False))
        .assign(batch=lambda d: d.samplename.str.extract("_([123])", expand=False))
    )
    sample_df.index = sample_df.samplename

    dds = py_DESeq2(
        count_matrix=df,
        design_matrix=sample_df,
        design_formula="~ batch + sample",
        gene_column="id",
    )

    dds.run_deseq()
    dds.get_deseq_result()
    res = dds.deseq_result
    assert res.query("padj < 0.05").shape == (35, 7)

    dds.get_deseq_result(contrast=["sample", "B", "A"])
    res = dds.deseq_result
    assert res.query("padj < 0.05").shape == (35, 7)

    lfc_res = dds.lfcShrink(coef=4, method="apeglm")
    assert lfc_res[lfc_res.padj < 0.05].shape[0] == 35

    res.to_csv(test_data_path + "/py_deseq.tsv", index=False, sep="\t")

    dds.run_deseq(test="LRT", reduced="~ batch")
    dds.get_deseq_result()
    res = dds.deseq_result
    res.to_csv(test_data_path + "/py_deseq_reduced.tsv", index=False, sep="\t")


@pytest.mark.parametrize(
    "case,r_table,py_table",
    [
        ("deseq", "R_deseq.tsv", "py_deseq.tsv"),
        ("deseq reduced", "R_deseq_reduced.tsv", "py_deseq_reduced.tsv"),
    ],
)
def test_result(run_r, case, r_table, py_table):
    os.chdir(os.path.dirname(test_data_path))

    py_tab = pd.read_table(os.path.join(test_data_path, py_table))
    r_tab = pd.read_table(os.path.join(test_data_path, r_table))

    for col in py_tab.columns:
        if py_tab.columns.dtype == "float64":
            assert np.all(
                np.isclose(py_tab[col].fillna(0), r_tab[col].fillna(0))
            ), f"{case} failed"
