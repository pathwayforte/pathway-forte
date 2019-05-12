# -*- coding: utf-8 -*-

"""CLI wrapper to perform GSEA with MSigDB gene sets."""

import logging
import os

import pandas as pd

from pathway_forte.constants import (
    CONCATENATED_MERGE_GENE_SETS, CONCATENATED_MERGE_GSEA_TSV, EXPRESSION_MATRIX, KEGG_MSIG_GSEA_TSV,
    KEGG_MSIG_SSGSEA_TSV, MERGE_GSEA, MERGE_SSGSEA, MSIGDB_KEGG_GENE_SETS, MSIGDB_REACTOME_GENE_SETS, MSIG_GSEA,
    MSIG_SSGSEA, NORMAL_EXPRESSION_SAMPLES, PHENOTYPE_CLASSES, REACTOME_MSIG_GSEA_TSV, REACTOME_MSIG_SSGSEA_TSV, TODAY,
    TUMOR_EXPRESSION_SAMPLES, make_gsea_export_directories, make_ssgsea_export_directories,
)
from pathway_forte.pathway_enrichment.functional_class import (
    create_cls_file, filter_gene_exp_data, run_gsea, run_ssgsea,
)

__all__ = [
    'do_gsea_msig',
]

logger = logging.getLogger(__name__)


def do_gsea_msig(data):
    make_gsea_export_directories()
    make_ssgsea_export_directories()

    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')
    num_normal_samples = NORMAL_EXPRESSION_SAMPLES.format(data)
    num_tumor_samples = TUMOR_EXPRESSION_SAMPLES.format(data)

    phenotype_path = PHENOTYPE_CLASSES.format(data)

    if not os.path.isfile(phenotype_path):
        logger.info('Creating cls file')
        create_cls_file(gene_exp, num_normal_samples, num_tumor_samples, data)

    logger.info('Running KEGG')
    kegg_gsea_results = run_gsea(
        gene_exp, MSIGDB_KEGG_GENE_SETS, phenotype_path, permutations=100, output_dir=MSIG_GSEA
    )
    kegg_gsea_results.res2d.to_csv(KEGG_MSIG_GSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running Reactome')
    reactome_gsea_results = run_gsea(
        gene_exp, MSIGDB_REACTOME_GENE_SETS, phenotype_path, permutations=100, output_dir=MSIG_GSEA
    )
    reactome_gsea_results.res2d.to_csv(REACTOME_MSIG_GSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running Concatenated Merge')
    conca_merge_gsea_results = run_gsea(
        gene_exp, CONCATENATED_MERGE_GENE_SETS, phenotype_path, permutations=100, output_dir=MERGE_GSEA
    )
    conca_merge_gsea_results.res2d.to_csv(CONCATENATED_MERGE_GSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running KEGG')
    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, MSIGDB_KEGG_GENE_SETS)
    kegg_ssgsea_results = run_ssgsea(filtered_expression_data, MSIGDB_KEGG_GENE_SETS, output_dir=MSIG_SSGSEA)
    kegg_ssgsea_results.res2d.to_csv(KEGG_MSIG_SSGSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running Reactome')
    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, MSIGDB_REACTOME_GENE_SETS)
    reactome_ssgsea_results = run_ssgsea(filtered_expression_data, MSIGDB_REACTOME_GENE_SETS, output_dir=MSIG_SSGSEA)
    reactome_ssgsea_results.res2d.to_csv(REACTOME_MSIG_SSGSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running Concatenated Merge')
    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, CONCATENATED_MERGE_GENE_SETS)
    conca_merge_ssgsea_results = run_ssgsea(
        filtered_expression_data, CONCATENATED_MERGE_GENE_SETS, output_dir=MERGE_SSGSEA
    )
    conca_merge_ssgsea_results.res2d.to_csv(CONCATENATED_MERGE_GSEA_TSV.format(data, TODAY), sep='\t')
