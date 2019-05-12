# -*- coding: utf-8 -*-

"""CLI wrapper to perform binary classification."""

import logging
import os

import pandas as pd

from pathway_forte.constants import (
    CLASSIFIER_RESULTS, EXPRESSION_MATRIX, KEGG, KEGG_SSGSEA, MERGE_SSGSEA, NORMAL_EXPRESSION_SAMPLES,
    PHENOTYPE_CLASSES, REACTOME, REACTOME_SSGSEA, TUMOR_EXPRESSION_SAMPLES, WIKIPATHWAYS, WIKIPATHWAYS_SSGSEA,
    make_classifier_results_directory,
)
from pathway_forte.pathway_enrichment.functional_class import create_cls_file
from pathway_forte.prediction.binary import get_l1_ratios, ssgsea_nes_to_df, train_elastic_net_model
from pathway_forte.utils import plot_aucs

__all__ = [
    'do_binary_prediction',
]

logger = logging.getLogger(__name__)


def do_binary_prediction(data, outer_cv_splits, inner_cv_splits, max_iter):
    make_classifier_results_directory()
    logger.info(f'Training Elastic Net via nested CV in the {data} dataset with'
                f' {outer_cv_splits} outer loops and {inner_cv_splits} inner loops')

    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')
    num_normal_samples = NORMAL_EXPRESSION_SAMPLES.format(data)
    num_tumor_samples = TUMOR_EXPRESSION_SAMPLES.format(data)

    if not os.path.isfile(PHENOTYPE_CLASSES.format(data)):
        logger.info('Creating cls file')
        create_cls_file(gene_exp, num_normal_samples, num_tumor_samples, data)

    phenotypes = PHENOTYPE_CLASSES.format(data)

    l1_ratio = get_l1_ratios()
    logger.info(f'L1 ratios: {l1_ratio}')

    kegg_ssgsea_nes_path = os.path.join(KEGG_SSGSEA, f'kegg_{data}.tsv')
    reactome_ssgsea_nes_path = os.path.join(REACTOME_SSGSEA, f'reactome_{data}.tsv')
    wikipathways_ssgsea_nes_path = os.path.join(WIKIPATHWAYS_SSGSEA, f'wikipathways_{data}.tsv')
    merge_ssgsea_nes_path = os.path.join(MERGE_SSGSEA, f'merge_{data}.tsv')
    _ds = [
        ('KEGG', kegg_ssgsea_nes_path, KEGG),
        ('Reactome', reactome_ssgsea_nes_path, REACTOME),
        ('WikiPathways', wikipathways_ssgsea_nes_path, WIKIPATHWAYS),
        ('Merge', merge_ssgsea_nes_path, 'MERGE')
    ]
    for name, tsv_path, fmt in _ds:
        logger.info(f'Training on {name}')
        x, y = ssgsea_nes_to_df(tsv_path, phenotypes)
        aucs = train_elastic_net_model(
            x,
            y,
            outer_cv_splits=outer_cv_splits,
            inner_cv_splits=inner_cv_splits,
            l1_ratio=l1_ratio,
            model_name=f"{data}-{fmt}",
            max_iter=max_iter,
        )
        plot_aucs(aucs, data, fmt, CLASSIFIER_RESULTS)
        logger.info(f"{name} AUCs: {aucs}")
