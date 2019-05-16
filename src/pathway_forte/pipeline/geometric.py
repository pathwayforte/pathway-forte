# -*- coding: utf-8 -*-

"""CLI wrapper to perform ORA using one-tailed hyper-geometric tests."""

import json
import logging
import os

from gseapy.parser import gsea_gmt_parser

from pathway_forte.pathway_enrichment.over_representation import (
    filter_p_value, perform_hypergeometric_test, read_fold_change_df,
)

__all__ = [
    'do_geometric',
]

logger = logging.getLogger(__name__)


def do_geometric(genesets, fold_changes, threshold):
    fc_df = read_fold_change_df(fold_changes)

    significant_genes = filter_p_value(fc_df, threshold)

    logger.info(f'There are a total of {len(significant_genes)} significant genes from the original {fc_df.shape[0]}.')

    # Note that the parser filters out gene sets smaller than 3 and larger than 1000
    gene_sets = gsea_gmt_parser(genesets)

    enriched_pathways = perform_hypergeometric_test(
        significant_genes,
        gene_sets,
        apply_threshold=threshold,
    )

    if threshold:
        logger.info('Filtering out pathways with q values > 0.05 according to fdr_bh')

    output = os.path.join(os.getcwd(), 'results.json')
    # Export dictionary as JSON
    with open(output, 'w') as file:
        json.dump(enriched_pathways, file, sort_keys=True, indent=2)

    logger.info(f'Results exported to {output}. # of pathways enriched {len(enriched_pathways)}')
