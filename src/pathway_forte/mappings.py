# -*- coding: utf-8 -*-

"""Function to deal with ComPath mappings."""

from collections import defaultdict
from functools import partial
from typing import Iterable, List, Mapping, Optional, Tuple

import numpy as np
import pandas as pd
from bio2bel.downloading import make_df_getter
from scipy.stats import ks_2samp, kstest, wilcoxon

from .constants import (
    IS_PART_OF, KEGG, KEGG_REACTOME_PATH, KEGG_REACTOME_URL, KEGG_WP_PATH, KEGG_WP_URL, MAPPING_TYPE, SOURCE_ID,
    SOURCE_RESOURCE, SPECIAL_MAPPINGS_PATH, SPECIAL_MAPPINGS_URL, TARGET_ID, TARGET_RESOURCE, WP_REACTOME_PATH,
    WP_REACTOME_URL,
)

__all__ = [
    'get_mapping_dict',
    'get_equivalent_pairs',
    'load_compath_mapping_dfs',
    'get_equivalent_mappings_dict',
    'remap_comparison_df',
    'get_equivalent_mapping_paired_test',
    'get_mlp_distribution_tests',
]

Identifier = Tuple[str, str]
EquivalenceMapping = Mapping[Identifier, List[Identifier]]

get_kegg_reactome_df = make_df_getter(KEGG_REACTOME_URL, KEGG_REACTOME_PATH)
get_wp_reactome_df = make_df_getter(WP_REACTOME_URL, WP_REACTOME_PATH)
get_kegg_wp_df = make_df_getter(KEGG_WP_URL, KEGG_WP_PATH)
get_special_mappings_df = make_df_getter(SPECIAL_MAPPINGS_URL, SPECIAL_MAPPINGS_PATH)


def get_mapping_dict(df: pd.DataFrame, mapping_type: str) -> Mapping[Identifier, List[Identifier]]:
    """Create a dictionary with ComPath mappings for each pathway."""
    mapping_dict = defaultdict(list)

    for index, row in df.iterrows():
        if row[MAPPING_TYPE] != mapping_type:
            continue

        if row[SOURCE_RESOURCE] != KEGG and row[TARGET_RESOURCE] != KEGG:
            mapping_dict[(row[SOURCE_RESOURCE], row[SOURCE_ID])].append((row[TARGET_RESOURCE], row[TARGET_ID]))
            mapping_dict[(row[TARGET_RESOURCE], row[TARGET_ID])].append((row[SOURCE_RESOURCE], row[SOURCE_ID]))

        elif row[SOURCE_RESOURCE] == KEGG and row[TARGET_RESOURCE] == KEGG:
            mapping_dict[(row[SOURCE_RESOURCE], row[SOURCE_ID].replace("path:", ""))].append(
                (row[TARGET_RESOURCE], row[TARGET_ID].replace("path:", "")))
            mapping_dict[(row[TARGET_RESOURCE], row[TARGET_ID].replace("path:", ""))].append(
                (row[SOURCE_RESOURCE], row[SOURCE_ID].replace("path:", "")))

        else:
            mapping_dict[(row[SOURCE_RESOURCE], row[SOURCE_ID].replace("path:", ""))].append(
                (row[TARGET_RESOURCE], row[TARGET_ID]))
            mapping_dict[(row[TARGET_RESOURCE], row[TARGET_ID])].append(
                (row[SOURCE_RESOURCE], row[SOURCE_ID].replace("path:", "")))

    return dict(mapping_dict)


def get_equivalent_pairs(df: pd.DataFrame):
    """Get equivalent pairs of pathways from 2 databases.

    :param df: pairwise mappings dataframe
    :return: equivalent pathway pairs dictionary {(SOURCE_RESOURCE,SOURCE_ID):[(TARGET_RESOURCE,TARGET_ID)]}
    :rtype: dict[list]
    """
    equivalent_pairs_dict = defaultdict(list)

    for index, row in df.iterrows():
        if row[MAPPING_TYPE] == IS_PART_OF:
            continue

        if row[SOURCE_RESOURCE] != KEGG and row[TARGET_RESOURCE] != KEGG:
            equivalent_pairs_dict[(row[SOURCE_RESOURCE], row[SOURCE_ID])].append(
                (row[TARGET_RESOURCE], row[TARGET_ID])
            )

        elif row[SOURCE_RESOURCE] == KEGG and row[TARGET_RESOURCE] == KEGG:
            equivalent_pairs_dict[(row[SOURCE_RESOURCE], row[SOURCE_ID].replace("path:", ""))].append(
                (row[TARGET_RESOURCE], row[TARGET_ID].replace("path:", ""))
            )

        else:
            equivalent_pairs_dict[(row[SOURCE_RESOURCE], row[SOURCE_ID].replace("path:", ""))].append(
                (row[TARGET_RESOURCE], row[TARGET_ID])
            )

    return dict(equivalent_pairs_dict)


def load_compath_mapping_dfs() -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load ComPath mappings data frames."""
    kegg_reactome_df = get_kegg_reactome_df()
    kegg_wikipathways_df = get_kegg_wp_df()
    wikipathways_reactome_df = get_wp_reactome_df()
    special_mappings_df = get_special_mappings_df()

    return (
        kegg_reactome_df,
        kegg_wikipathways_df,
        wikipathways_reactome_df,
        special_mappings_df,
    )


def get_equivalent_mappings_dict() -> EquivalenceMapping:
    """Get mapping dictionary of all equivalent pairs of pathways.

    Special mappings are not included in the overall mappings as some of the WP pathways possess identical IDs.
    """
    kegg_reactome_df, kegg_wikipathways_df, wikipathways_reactome_df, special_mappings_df = load_compath_mapping_dfs()

    # Get mapping dictionary of all equivalent pairs of pathways. Special mappings are not included in the overall
    # mappings as some of the WP pathways possess identical IDs.
    return get_mapping_dict(
        pd.concat([kegg_reactome_df, kegg_wikipathways_df, wikipathways_reactome_df]),
        'equivalentTo'
    )


def get_wikipathways(pathways):
    return {
        pathway: pathway
        for pathway in pathways
        if str(pathway).startswith('WP')
    }


def get_reactomes(pathways):
    return {
        pathway: f"R-HSA-{pathway}"
        for pathway in pathways
        if not str(pathway).startswith('WP') and len(str(pathway)) > 4
    }


def get_keggs(pathways):
    return {
        pathway: f"hsa0{pathway}"
        for pathway in pathways
        if not str(pathway).startswith('WP') and len(str(pathway)) <= 4
    }


def get_equivalent_mapping_paired_test(
        df: pd.DataFrame,
        source_db: str,
        target_db: str,
        *,
        db_column_name: str,
        identifier_column_name: str,
        pval_column_name: str,
        equivalent_mappings_dict: Optional[EquivalenceMapping] = None,
        test: str = 'wilcoxon',
) -> float:
    remapped_df = remap_comparison_df(
        df=df,
        source_db=source_db,
        target_db=target_db,
        db_column_name=db_column_name,
        identifier_column_name=identifier_column_name,
        pval_column_name=pval_column_name,
        equivalent_mappings_dict=equivalent_mappings_dict,
    )

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
        db_column_name: str,
        identifier_column_name: str,
        pval_column_name: str,
        equivalent_mappings_dict: Optional[EquivalenceMapping] = None,
) -> pd.DataFrame:
    """
    
    :param df: 
    :param source_db: This is the database that becomes the left one
    :param target_db:
    :param equivalent_mappings_dict: 
    """
    if equivalent_mappings_dict is None:
        equivalent_mappings_dict = get_equivalent_mappings_dict()

    source_indexes = df[db_column_name] == source_db

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
        db_column_name: str,
        equivalent_mappings_dict: Optional[EquivalenceMapping] = None,
        alpha: float = 0.05,
) -> pd.DataFrame:
    rows = []

    for dataset in datasets:
        dataset_df = df[df.dataset == dataset]
        for comparison in dataset_df.comparison.unique():
            dataset_comparison_df = dataset_df[dataset_df.comparison == comparison]

            # 2 sample KS test
            _, ks_p = ks_2samp(*[
                group[pval_column_name].values
                for name, group in dataset_comparison_df.groupby(db_column_name)
            ])

            # Wilcoxon test of the differences of the paired differences of MLP-values
            source_db, target_db = comparison.split('_')
            wilcoxon_p = get_equivalent_mapping_paired_test(
                dataset_comparison_df,
                source_db,
                target_db,
                identifier_column_name=identifier_column_name,
                pval_column_name=pval_column_name,
                db_column_name=db_column_name,
                equivalent_mappings_dict=equivalent_mappings_dict,
                test='wilcoxon',
            )

            # KS Test of the differences of the paired differences of MLP-values
            ks_paired_mlp_p = get_equivalent_mapping_paired_test(
                dataset_comparison_df,
                source_db,
                target_db,
                identifier_column_name=identifier_column_name,
                pval_column_name=pval_column_name,
                db_column_name=db_column_name,
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
