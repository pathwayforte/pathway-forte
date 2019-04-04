# -*- coding: utf-8 -*-

"""This module contains the functions to run Over Representation Analysis (ORA)."""
import itertools as itt
import logging

import gseapy
import numpy as np

log = logging.getLogger(__name__)


def _prepare_hypergeometric_test(query_gene_set, pathway_gene_set, gene_universe):
    """Prepare the matrix for hypergeometric test calculations.

    :param set[str] query_gene_set: gene set to test against pathway
    :param set[str] pathway_gene_set: pathway gene set
    :param int gene_universe: number of HGNC symbols
    :rtype: numpy.ndarray
    :return: 2x2 matrix
    """
    return np.array(
        [[len(query_gene_set.intersection(pathway_gene_set)),
          len(query_gene_set.difference(pathway_gene_set))
          ],
         [len(pathway_gene_set.difference(query_gene_set)),
          gene_universe - len(pathway_gene_set.union(query_gene_set))
          ]
         ]
    )


def filter_gene_exp_data(expression_data, gmt_file):
    """Filter gene expression data file to include only gene names which are found in the gene set files.

    :param pandas.core.frame.DataFrame expression_data: gene expression values for samples
    :param str gmt_file: .gmt file containing gene sets
    :return: Filtered gene expression data with genes with no correspondences in gene sets removed
    :rtype: pandas.core.frame.DataFramekegg_xml_parser.py
    """
    filtered_expression_data = expression_data.copy()

    # Gene universe from gene set
    gene_sets = gseapy.parser.gsea_gmt_parser(gmt_file, max_size=40000)

    gene_universe = set(itt.chain(*gene_sets.values()))

    counter = 0

    genes_to_remove = list()

    for gene in filtered_expression_data.index.values:

        if gene not in gene_universe:
            genes_to_remove.append(gene)
            counter += 1

    log.info(f'Expression data has {len(filtered_expression_data.index.values)}')
    log.info(f'Gene universe has {len(gene_universe)}')
    log.info(f'{counter} were removed in expression data')
    log.info(
        f'{(len(filtered_expression_data.index.values)-counter) /len(gene_universe)*100:.4f}% '
        f'of the gene expression data is mapped to the pathway datasets')

    return filtered_expression_data.drop(genes_to_remove)
