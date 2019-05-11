# -*- coding: utf-8 -*-

"""Function to deal with ComPath mappings."""

from collections import defaultdict
from typing import Tuple

import pandas as pd

from pathway_forte.constants import (
    IS_PART_OF, KEGG, KEGG_REACTOME_URL, KEGG_WP_URL, MAPPING_TYPE, SOURCE_ID,
    SOURCE_RESOURCE, SPECIAL_MAPPINGS_URL, TARGET_ID, TARGET_RESOURCE, WP_REACTOME_URL,
)

__all__ = [
    'get_mapping_dict',
    'get_equivalent_pairs',
    'load_compath_mapping_dfs',
    'get_equivalent_mappings_dict',
]


def get_mapping_dict(df: pd.DataFrame, mapping_type: str):
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

    return mapping_dict


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
    """Load ComPath mappings dataframes."""
    kegg_reactome_df = pd.read_csv(KEGG_REACTOME_URL)
    kegg_wikipathways_df = pd.read_csv(KEGG_WP_URL)
    wikipathways_reactome_df = pd.read_csv(WP_REACTOME_URL)
    special_mappings_df = pd.read_csv(SPECIAL_MAPPINGS_URL)
    return (
        kegg_reactome_df,
        kegg_wikipathways_df,
        wikipathways_reactome_df,
        special_mappings_df,
    )


def get_equivalent_mappings_dict():
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
