# -*- coding: utf-8 -*-

"""CLI wrapper to perform GSEA with the four major pathway databases."""

import logging
import os

import pandas as pd

from pathway_forte.constants import (
    EXPRESSION_MATRIX, KEGG_GSEA, KEGG_GSEA_TSV, MERGE_GSEA, MERGE_GSEA_TSV, NORMAL_EXPRESSION_SAMPLES,
    PHENOTYPE_CLASSES, REACTOME_GSEA, REACTOME_GSEA_TSV, TODAY, TUMOR_EXPRESSION_SAMPLES, WIKIPATHWAYS_GSEA,
    WIKIPATHWAYS_GSEA_TSV, check_gmt_files, make_gsea_export_directories,
)
from pathway_forte.pathway_enrichment.functional_class import create_cls_file, run_gsea

__all__ = [
    'do_gsea',
]

logger = logging.getLogger(__name__)


def do_gsea(data, permutations):
    make_gsea_export_directories()

    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')
    num_normal_samples = NORMAL_EXPRESSION_SAMPLES.format(data)
    num_tumor_samples = TUMOR_EXPRESSION_SAMPLES.format(data)

    phenotype_path = PHENOTYPE_CLASSES.format(data)

    if not os.path.isfile(phenotype_path):
        logger.info('Creating cls file')
        create_cls_file(gene_exp, num_normal_samples, num_tumor_samples, data)

    kegg_gene_set, reactome_gene_set, wikipathways_gene_set, _ = check_gmt_files()
    _ds = [
        ('KEGG', kegg_gene_set, KEGG_GSEA, KEGG_GSEA_TSV),
        ('Reactome', reactome_gene_set, REACTOME_GSEA, REACTOME_GSEA_TSV),
        ('WikiPathways', wikipathways_gene_set, WIKIPATHWAYS_GSEA, WIKIPATHWAYS_GSEA_TSV),
    ]
    for name, gene_set, output_dir, tsv_fmt in _ds:
        logger.info(f'Running {name}')
        results = run_gsea(
            gene_exp,
            gene_set,
            phenotype_path,
            permutations=permutations,
            output_dir=output_dir,
        )
        results.res2d.to_csv(tsv_fmt.format(data, TODAY), sep='\t')
