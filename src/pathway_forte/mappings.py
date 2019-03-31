# -*- coding: utf-8 -*-

"""Function to deal with ComPath mappings."""

from collections import defaultdict

import pandas as pd

from pathway_forte.constants import *


def get_mapping_dict(df, mapping_type):
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


def get_equivalent_pairs(df):
    """Get equivalent pairs of pathways from 2 databases.

    :param pandas.core.frame.DataFrame df: pairwise mappings dataframe
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

    return equivalent_pairs_dict


def load_compath_mapping_dfs():
    """Load ComPath mappings dfs"""
    kegg_reactome_df = pd.read_csv(KEGG_REACTOME_URL)
    kegg_wikipathways_df = pd.read_csv(KEGG_WP_URL)
    wikipathways_reactome_df = pd.read_csv(WP_REACTOME_URL)
    special_mappings_df = pd.read_csv(SPECIAL_MAPPINGS_URL)

    return kegg_reactome_df, kegg_wikipathways_df, wikipathways_reactome_df, special_mappings_df
