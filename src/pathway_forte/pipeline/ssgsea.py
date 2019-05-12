# -*- coding: utf-8 -*-

"""CLI wrapper to perform ssGSEA."""

import logging

import pandas as pd

from pathway_forte.constants import (
    EXPRESSION_MATRIX, KEGG_SSGSEA, KEGG_SSGSEA_TSV, MERGE_SSGSEA, MERGE_SSGSEA_TSV, REACTOME_SSGSEA,
    REACTOME_SSGSEA_TSV, TODAY, WIKIPATHWAYS_SSGSEA, WIKIPATHWAYS_SSGSEA_TSV, check_gmt_files,
    make_ssgsea_export_directories,
)
from pathway_forte.pathway_enrichment.functional_class import filter_gene_exp_data, run_ssgsea

__all__ = [
    'do_ssgsea',
]

logger = logging.getLogger(__name__)


def do_ssgsea(data):
    make_ssgsea_export_directories()

    # Read data
    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')

    kegg_gene_set, reactome_gene_set, wikipathways_gene_set, merge_gene_set = check_gmt_files()

    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, merge_gene_set)

    _ds = [
        ('KEGG', kegg_gene_set, KEGG_SSGSEA, KEGG_SSGSEA_TSV),
        ('Reactome', reactome_gene_set, REACTOME_SSGSEA, REACTOME_SSGSEA_TSV),
        ('WikiPathways', wikipathways_gene_set, WIKIPATHWAYS_SSGSEA, WIKIPATHWAYS_SSGSEA_TSV),
        ('MergeDataset', merge_gene_set, MERGE_SSGSEA, MERGE_SSGSEA_TSV),
    ]
    for name, gene_set, output_dir, fmt in _ds:
        logger.info(f'Running {name}')
        results = run_ssgsea(filtered_expression_data, gene_set, output_dir=output_dir)
        results.res2d.to_csv(fmt.format(data, TODAY), sep='\t')
