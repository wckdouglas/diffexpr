#!/usr/bin/env python
import os
import re
import shlex
import subprocess
import warnings

import numpy as np
import pandas as pd
import pytest

from diffexpr.py_deseq import py_DESeq2

warnings.filterwarnings("ignore")
test_data_path = os.path.dirname(os.path.realpath(__file__)) + "/data"


@pytest.fixture(scope="module", name="run_r")
def fixture_run_r():
    os.chdir(os.path.dirname(test_data_path))
    os.system("Rscript deseq.R")


def count_matrix():
    """
        id     A_1     A_2     A_3     B_1     B_2     B_3
    0  ERCC-00002  111461  106261  107547  333944  199252  186947
    1  ERCC-00003    6735    5387    5265   13937    8584    8596
    2  ERCC-00004   17673   13983   15462    5065    3222    3353
    3  ERCC-00009    4669    4431    4211    6939    4155    3647
    4  ERCC-00012       0       2       0       0       0       0
    """
    return pd.read_csv(test_data_path + "/ercc.tsv", sep="\t")


def sample_metadata():
    df = count_matrix()
    sample_df = (
        pd.DataFrame({"samplename": df.columns})
        .query('samplename != "id"')
        .assign(sample=lambda d: d["samplename"].str.extract("([AB])_", expand=False))
        .assign(batch=lambda d: d["samplename"].str.extract("_([123])", expand=False))
        .pipe(lambda d: d.set_index(d["samplename"]))
        # this duplicates the "samplename" column into index, needed for `setup_deseq`
    )
    return sample_df


@pytest.fixture(scope="module")
def setup_deseq():
    df = count_matrix()
    sample_df = sample_metadata()

    dds = py_DESeq2(
        count_matrix=df,
        design_matrix=sample_df,
        design_formula="~ batch + sample",
        gene_column="id",
    )

    return df, dds


@pytest.mark.parametrize("matrix", [("count"), ("metadata")])
def test_deseq_not_dataframe_exception(matrix):
    df = pd.read_csv(test_data_path + "/ercc.tsv", sep="\t")
    sample_df = sample_metadata()
    if matrix == "count":
        df = df.values
    elif matrix == "metadata":
        sample_df = sample_df.values

    with pytest.raises(ValueError, match="should be pd.DataFrame type"):
        py_DESeq2(
            count_matrix=df,
            design_matrix=sample_df,
            design_formula="~ batch + sample",
            gene_column="id",
        )


def test_deseq_no_id_exception():
    df = pd.read_csv(test_data_path + "/ercc.tsv", sep="\t")
    sample_df = sample_metadata()
    with pytest.raises(ValueError, match="The given gene_column name is not a column"):
        py_DESeq2(
            count_matrix=df,
            design_matrix=sample_df,
            design_formula="~ batch + sample",
            gene_column="id1",
        )


def test_deseq(setup_deseq):
    df, dds = setup_deseq
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


def test_normalized_count(setup_deseq):
    df, dds = setup_deseq
    norm_count_df = dds.normalized_count()
    assert norm_count_df.shape == df.shape


@pytest.mark.parametrize(
    "blind,fit_type",
    [(True, "parametric"), (True, "local"), (True, "mean"), (False, "parametric"), (False, "local"), (False, "mean")],
)
def test_vst(setup_deseq, blind, fit_type):
    df, dds = setup_deseq
    vst = dds.vst(blind=blind, fit_type=fit_type)
    assert vst.shape == df.shape


@pytest.mark.parametrize(
    "blind,fit_type",
    [(True, "parametric"), (True, "local"), (True, "mean"), (False, "parametric"), (False, "local"), (False, "mean")],
)
def test_rlog(setup_deseq, blind, fit_type):
    df, dds = setup_deseq
    rlog = dds.rlog(blind=blind, fit_type=fit_type)
    assert rlog.shape == df.shape


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
            assert np.all(np.isclose(py_tab[col].fillna(0), r_tab[col].fillna(0))), f"{case} failed"


def test_deseq2_version(setup_deseq):
    _, dds = setup_deseq

    # extract R deseq2 version
    re_ver = re.compile("\d+.\d+.\d+")
    cmd = """
        Rscript -e "packageVersion('DESeq2')" 
    """
    version_output = subprocess.check_output(shlex.split(cmd)).decode()
    version = re_ver.findall(version_output)[0]

    # let's see if we extracted the correct version
    assert dds.deseq2_version == version
