# -*- coding: utf-8 -*-

"""CLI wrapper to perform ORA using one-tailed hyper-geometric tests."""

import logging
import os
from typing import Optional

from gseapy.parser import gsea_gmt_parser

from pathway_forte.pathway_enrichment.over_representation import (
    filter_p_value, perform_hypergeometric_test, read_fold_change_df,
)

__all__ = [
    'do_geometric',
]

logger = logging.getLogger(__name__)


def do_geometric(
        gmt_path: str,
        fold_changes_path: str,
        apply_threshold: bool,
        output: Optional[str] = None,
):
    """Wrapper to run hyper-geometric test."""
    fc_df = read_fold_change_df(fold_changes_path)

    significant_genes = filter_p_value(fc_df, apply_threshold)

    logger.info(f'There are a total of {len(significant_genes)} significant genes from the original {fc_df.shape[0]}.')

    # Note that the parser filters out gene sets smaller than 3 and larger than 1000
    gene_sets = gsea_gmt_parser(gmt_path)

    df = perform_hypergeometric_test(
        significant_genes,
        gene_sets,
        apply_threshold=apply_threshold,
    )

    if apply_threshold:
        logger.info('Filtering out pathways with q values > 0.05 according to fdr_bh')

    if output is None:
        output = os.path.join(os.getcwd(), 'results.tsv')

    df.to_csv(output, sep='\t', index=False)
    logger.info(f'Results exported to {output}. # of pathways enriched {len(df.index)}')
