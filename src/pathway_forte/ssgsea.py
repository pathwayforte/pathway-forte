# -*- coding: utf-8 -*-

"""This module runs ssGSEA on the TCGA dataset."""

import itertools as itt
import logging

import gseapy

from .constants import SSGSEA

log = logging.getLogger(__name__)


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


def run_ssgsea(filtered_expression_data, gene_set, output_dir=SSGSEA, processes=1):
    """Run single sample GSEA (ssGSEA) on filtered gene expression data set.

    :param pandas.core.frame.DataFrame expression_data: filtered gene expression values for samples
    :param str gmt_file: .gmt file containing gene sets
    :param output_dir: output directory
    :return: ssGSEA results in respective directory
    """
    ssgsea_result = gseapy.ssgsea(
        data=filtered_expression_data,
        gene_sets=gene_set,
        outdir=output_dir,  # do not write output to disk
        max_size=3000,
        sample_norm_method='rank',  # choose 'custom' for your own rank list
        permutation_num=0,  # skip permutation procedure, because you don't need it
        no_plot=True,  # skip plotting to speed up
        processes=processes, format='png'
    )
    log.info('Done with ssGSEA')
    return ssgsea_result
