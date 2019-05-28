# -*- coding: utf-8 -*-

"""This module contains the functions to run Over Representation Analysis (ORA)."""

import logging
from typing import Mapping, Set

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from pathway_forte.constants import FC_COLUMNS, GENE_SYMBOL, P_VALUE

logger = logging.getLogger(__name__)


def read_fold_change_df(path: str) -> pd.DataFrame:
    """Read csv with gene names, fold changes and their p-values."""
    df = pd.read_csv(path)

    # Check all columns are present
    missing_columns = [
        column
        for column in FC_COLUMNS
        if column not in df
    ]
    if missing_columns:
        raise ValueError(f'Missing columns {", ".join(missing_columns)} in {path}')

    return df


def filter_p_value(df: pd.DataFrame, p_value: bool = True, cutoff: float = 0.05):
    """Return significantly differentially expressed genes in fold change df."""
    # Apply p_value threshold too
    if p_value:
        query_exp = f'{P_VALUE} <= {cutoff}'
        df = df.query(query_exp)

    return set(df[GENE_SYMBOL])


def _prepare_hypergeometric_test(
        query_gene_set: Set[str],
        pathway_gene_set: Set[str],
        gene_universe: int,
) -> np.ndarray:
    """Prepare the matrix for hypergeometric test calculations.

    :param query_gene_set: gene set to test against pathway
    :param pathway_gene_set: pathway gene set
    :param gene_universe: number of HGNC symbols
    :return: 2x2 matrix
    """
    # Cast lists to sets
    if not isinstance(query_gene_set, set):
        query_gene_set = set(query_gene_set)
    if not isinstance(pathway_gene_set, set):
        pathway_gene_set = set(pathway_gene_set)

    # Return matrix to test hyper-geometric test
    return np.array([
        [
            len(query_gene_set.intersection(pathway_gene_set)),
            len(query_gene_set.difference(pathway_gene_set)),
        ],
        [
            len(pathway_gene_set.difference(query_gene_set)),
            gene_universe - len(pathway_gene_set.union(query_gene_set)),
        ],
    ])


def perform_hypergeometric_test(
        genes_to_test: Set[str],
        pathway_dict: Mapping[str, Set[str]],
        gene_universe: int = 41714,
        apply_threshold: bool = False,
        threshold: float = 0.05,
) -> pd.DataFrame:
    """Perform hypergeometric tests.

    :param genes_to_test: gene set to test against pathway
    :param pathway_dict: pathway name to gene set
    :param gene_universe: number of HGNC symbols
    :param apply_threshold: return only significant pathways
    :param threshold: significance threshold (by default 0.05)
    """
    rows = []
    for (pathway_id, database), pathway_gene_set in pathway_dict.items():
        # Prepare the test table to conduct the fisher test
        test_table = _prepare_hypergeometric_test(genes_to_test, pathway_gene_set, gene_universe)
        # Calculate fisher test (returns tuple of odds ratio and p_value
        p_value = fisher_exact(test_table, alternative='greater')[1]
        rows.append((database, pathway_id, p_value))

    df = pd.DataFrame(rows, columns=['database', 'pathway_id', 'pval'])
    correction_test = multipletests(df.pval, method='fdr_bh')
    df['qval'] = correction_test[1]

    if apply_threshold:
        logger.info('Filtering out pathways with q-values > 0.05 according to fdr_bh')
        df = df[df['qval'] < threshold]

    return df
