# -*- coding: utf-8 -*-

"""CLI wrapper to perform ssGSEA."""

import logging

from pathway_forte.constants import make_classifier_results_directory
from pathway_forte.prediction.binary import get_l1_ratios, train_elastic_net_model

__all__ = [
    'do_stability_prediction',
]

logger = logging.getLogger(__name__)


def do_stability_prediction(
        ssgsea_scores_path,
        phenotypes_path,
        *,
        outer_cv_splits,
        inner_cv_splits,
        max_iter,
):
    """Train elastic net."""
    make_classifier_results_directory()
    logger.info(
        f'Training Elastic Net via nested CV for {outer_cv_splits} dataset '
        f'with {ssgsea_scores_path} (phenotypes: {phenotypes_path}) '
        f'outer loops and {inner_cv_splits} inner loops'
    )

    l1_ratio = get_l1_ratios()
    logger.info(f'L1 ratios: {l1_ratio}')

    results = train_elastic_net_model(
        ssgsea_scores_path,
        phenotypes_path,
        outer_cv_splits=outer_cv_splits,
        inner_cv_splits=inner_cv_splits,
        l1_ratio=l1_ratio,
        model_name='elastic_net',
        max_iter=max_iter,
    )
    logger.info(results)
