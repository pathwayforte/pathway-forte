# -*- coding: utf-8 -*-

from functools import partial
from typing import Iterable, Optional

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, kstest, wilcoxon

from pathway_forte.mappings import EquivalenceMapping, get_equivalent_mappings_dict

__all__ = [
    'get_equivalent_mapping_paired_test',
    'get_mlp_distribution_tests',
]


def get_equivalent_mapping_paired_test(
        df: pd.DataFrame,
        source_db: str,
        target_db: str,
        *,
        database_column_name: str,
        identifier_column_name: str,
        pval_column_name: str,
        equivalent_mappings_dict: Optional[EquivalenceMapping] = None,
        test: str = 'wilcoxon',
) -> float:
    remapped_df = remap_comparison_df(
        df=df,
        source_db=source_db,
        target_db=target_db,
        database_column_name=database_column_name,
        identifier_column_name=identifier_column_name,
        pval_column_name=pval_column_name,
        equivalent_mappings_dict=equivalent_mappings_dict,
    )
    if 0 == len(remapped_df):
        raise ValueError(f'No comparison for {source_db} and {target_db}')

    if test == 'wilcoxon':
        test = wilcoxon
    elif test == 'ks':
        test = partial(kstest, cdf='uniform')
    else:
        raise ValueError('invalid test type: {test}')

    _, p_value = test(remapped_df['mlp_diff'].tolist())
    return p_value


def remap_comparison_df(
        df: pd.DataFrame,
        source_db: str,
        target_db: str,
        *,
        database_column_name: str,
        identifier_column_name: str,
        pval_column_name: str,
        equivalent_mappings_dict: Optional[EquivalenceMapping] = None,
) -> pd.DataFrame:
    if equivalent_mappings_dict is None:
        equivalent_mappings_dict = get_equivalent_mappings_dict()

    source_indexes = df[database_column_name] == source_db

    rows = []
    it = df.loc[source_indexes, [identifier_column_name, pval_column_name]].iterrows()
    for i, (source_identifier, source_p) in it:
        target_identifier, target_p = None, 1.0
        for x, y in equivalent_mappings_dict[source_db, source_identifier]:
            if target_db == x:
                target_identifier = y

                pval_rows = df[df[identifier_column_name] == target_identifier]
                try:
                    target_p = pval_rows.iloc[0][pval_column_name]
                except IndexError:  # if no p-value found, default to 1.0
                    target_p = 1.0

        source_mlp = -np.log10(source_p)
        target_mlp = -np.log10(target_p)
        rows.append((
            source_db,
            source_identifier,
            source_p,
            source_mlp,
            target_db,
            target_identifier,
            target_p,
            target_mlp,
            target_mlp - source_mlp,
        ))

    return pd.DataFrame(
        rows,
        columns=[
            'source_db',
            'source_identifier',
            'source_p',
            'source_mlp',
            'target_db',
            'target_identifier',
            'target_p',
            'target_mlp',
            'mlp_diff',
        ]
    )


def get_mlp_distribution_tests(
        df: pd.DataFrame,
        datasets: Iterable[str],
        *,
        identifier_column_name: str,
        pval_column_name: str,
        database_column_name: str,
        dataset_column_name: str,
        equivalent_mappings_dict: Optional[EquivalenceMapping] = None,
        alpha: float = 0.05,
) -> pd.DataFrame:
    _df = df[df.comparison.notna()]

    rows = []

    for dataset in datasets:
        dataset_df = _df[_df[dataset_column_name] == dataset]
        for comparison in dataset_df.comparison.unique():
            dataset_comparison_df = dataset_df[dataset_df.comparison == comparison]

            # 2 sample KS test
            groups = [
                group[pval_column_name].values
                for name, group in dataset_comparison_df.groupby(database_column_name)
            ]

            try:
                a, b = groups
            except ValueError:
                raise ValueError(f'problem with {dataset} and {comparison}')

            if 0 == len(a) or 0 == len(b):
                raise ValueError(f'problem with {dataset} and {comparison}')

            _, ks_p = ks_2samp(a, b)

            # Wilcoxon test of the differences of the paired differences of MLP-values
            source_db, target_db = comparison.split('_')
            try:
                wilcoxon_p = get_equivalent_mapping_paired_test(
                    dataset_comparison_df,
                    source_db,
                    target_db,
                    identifier_column_name=identifier_column_name,
                    pval_column_name=pval_column_name,
                    database_column_name=database_column_name,
                    equivalent_mappings_dict=equivalent_mappings_dict,
                    test='wilcoxon',
                )
            except ValueError:
                print(f'no comparison for {source_db} and {target_db} for {dataset}')
                continue

            # KS Test of the differences of the paired differences of MLP-values
            ks_paired_mlp_p = get_equivalent_mapping_paired_test(
                dataset_comparison_df,
                source_db,
                target_db,
                identifier_column_name=identifier_column_name,
                pval_column_name=pval_column_name,
                database_column_name=database_column_name,
                equivalent_mappings_dict=equivalent_mappings_dict,
                test='ks',
            )

            rows.append((
                dataset,
                comparison,
                ks_p,
                -np.log10(ks_p),
                wilcoxon_p,
                -np.log10(wilcoxon_p),
                ks_paired_mlp_p,
                -np.log10(ks_paired_mlp_p),
            ))

    rv = pd.DataFrame(
        rows,
        columns=[
            'dataset',
            'comparison',
            'ks_p',
            'ks_mlp',
            'wilcoxon_p',
            'wilcoxon_mlp',
            'ks_paired_p',
            'ks_paired_mlp',
        ],
    )
    rv['wilcoxon_significant'] = rv['wilcoxon_p'] < alpha
    rv['ks_paired_significant'] = rv['ks_paired_p'] < alpha
    return rv
