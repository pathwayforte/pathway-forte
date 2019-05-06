# -*- coding: utf-8 -*-

"""Test ORA import."""

import os
import unittest

from gseapy.parser import gsea_gmt_parser

from pathway_forte.pathway_enrichment.over_representation import (
    read_fold_change_df,
    filter_fold_change_fd,
    perform_hypergeometric_test
)

TEST_FOLDER = os.path.dirname(os.path.realpath(__file__))
GMT_EXAMPLE = os.path.join(TEST_FOLDER, 'data', 'example.gmt')
FOLD_CHANGES_EXAMPLE = os.path.join(TEST_FOLDER, 'data', 'example_fold_changes.csv')


class TestOra(unittest.TestCase):

    def test_get_df(self):
        """Test getting the genes out of the csv file."""
        fc_df = read_fold_change_df(FOLD_CHANGES_EXAMPLE)

        significant_genes = filter_fold_change_fd(fc_df)
        gene_sets = gsea_gmt_parser(GMT_EXAMPLE)

        self.assertEqual(significant_genes, {'C', 'A'})
        self.assertEqual(gene_sets, {'pathway1': ['A', 'B', 'C', 'D'], 'pathway2': ['E', 'F', 'G', 'H']})

        enriched_pathways = perform_hypergeometric_test(
            significant_genes,
            gene_sets,
            apply_threshold=True
        )

        self.assertEqual(len(enriched_pathways.keys()), 1)
        self.assertIn('pathway1', enriched_pathways.keys())
