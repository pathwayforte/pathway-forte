# -*- coding: utf-8 -*-

"""Test merge gene set."""

import logging
import unittest
from collections import Counter

import pandas as pd

from pathway_forte.mappings import load_compath_mapping_dfs, get_mapping_dict, check_gmt_files

logger = logging.getLogger(__name__)

# Mappings linked with duplicate WikiPathways/KEGG representation need to be skipped since they one pathway maps to
# multiple mappings
BLACK_LIST = {
    'R-HSA-5683057', 'WP623', 'WP61', 'WP3845', 'WP366', 'R-HSA-9006936', 'R-HSA-2028269', 'WP75', 'WP3858',
    'WP2586', 'WP2873', 'R-HSA-174403', 'WP100', 'WP4506', 'R-HSA-196819', 'WP4297', 'R-HSA-1430728', 'R-HSA-163210',
    'R-HSA-1430728', 'hsa00190', 'hsa04010', 'WP382', 'WP422', 'hsa04350', 'WP560', 'hsa04390', 'hsa04392'
}


class TestMergeGmt(unittest.TestCase):

    def test_gmt_file(self):

        kegg_reactome_df, kegg_wikipathways_df, wikipathways_reactome_df, special_mappings_df = load_compath_mapping_dfs()

        equivalent_mappings_dict = get_mapping_dict(
            pd.concat([kegg_reactome_df, kegg_wikipathways_df, wikipathways_reactome_df]),
            'equivalentTo'
        )

        _, _, _, merge_gene_set = check_gmt_files()

        with open(merge_gene_set) as f:
            content = f.readlines()
            # Get the two first cells in each row (pathway ids and resources)
            pathway_tuples = [
                line.split('\t')[0:2]
                for line in content
            ]

        # Zip ids and resources to a common
        pathways_mapped = [
            list(zip(resources.split('|'), pathway_ids.split('|')))
            for pathway_ids, resources in pathway_tuples
        ]

        for pathways in pathways_mapped:

            if 1 == len(pathways):
                # If the pathway doesnt have a mapping in the GMT file but it is in the equivalent mapping dict
                # RAISE ERROR
                if pathways[0] in equivalent_mappings_dict:
                    raise ValueError(f'{pathways} should have a mapping')
                continue

            copy_pathway = pathways

            for index, (resource, pathway) in enumerate(pathways):
                # Skip multiduplicate pathways
                if pathway in BLACK_LIST:
                    continue

                mapping_pathways_in_iteration = copy_pathway.copy()
                mapping_pathways_in_iteration.pop(index)

                real_mappings = equivalent_mappings_dict[(resource, pathway)]

                self.assertEqual(set(mapping_pathways_in_iteration), set(real_mappings))

        counter = Counter([
            pathway
            for pathways in pathways_mapped
            for resource, pathway in pathways
        ])

        # All pathways should be only present once. Check for duplicates
        for pathway, count in counter.items():
            self.assertEqual(count, 1)
