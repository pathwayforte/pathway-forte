# -*- coding: utf-8 -*-

"""This module contains the functions to run Over Representation Analysis (ORA)."""

import logging
from typing import Set

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from pathway_forte.constants import FC_COLUMNS, GENE_SYMBOL

log = logging.getLogger(__name__)


def read_fold_change_df(path) -> pd.DataFrame:
    """Read csv with gene names, fold changes and their pvalues."""
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


def filter_p_value(df, p_value=True, cutoff=0.05):
    """Return significantly differentially expressed genes in fold change df."""
    # Apply p_value threshold too
    if p_value:

        query_exp = f'{P_VALUE} <= {cutoff}'

        df = df.query(query_exp)

    return set(df[GENE_SYMBOL])


def _prepare_hypergeometric_test(query_gene_set: Set[str], pathway_gene_set: Set[str],
                                 gene_universe: int) -> np.ndarray:
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
        pathway_dict,
        gene_universe=41714,
        apply_threshold=False,
        threshold=0.05,
):
    """Perform hypergeometric tests.

    :param genes_to_test: gene set to test against pathway
    :param dict[str,set] pathway_dict: pathway name to gene set
    :param int gene_universe: number of HGNC symbols
    :param Optional[bool] apply_threshold: return only significant pathways
    :param Optional[float] threshold: significance threshold (by default 0.05)
    :rtype: dict[str,dict[str,dict]]
    :return: manager_pathways_dict with p value info
    """
    pathway_to_p_value = dict()
    results = dict()

    for pathway_id, pathway_gene_set in pathway_dict.items():
        # Prepare the test table to conduct the fisher test
        test_table = _prepare_hypergeometric_test(genes_to_test, pathway_gene_set, gene_universe)

        # Calculate fisher test
        oddsratio, pvalue = fisher_exact(test_table, alternative='greater')

        pathway_to_p_value[pathway_id] = pvalue

    # Split the dictionary into names_id tuples and p values to keep the same order
    manager_pathway_id, p_values = zip(*pathway_to_p_value.items())
    correction_test = multipletests(p_values, method='fdr_bh')

    q_values = correction_test[1]

    # Update original dict with p value corrections
    for i, pathway_id in enumerate(manager_pathway_id):
        results[pathway_id] = q_value = round(q_values[i], 4)

        # [Optional] Delete the pathway if does not pass the threshold
        if apply_threshold and q_value > threshold:
            del results[pathway_id]

    return results
