# -*- coding: utf-8 -*-

""""""

import logging

import pandas as pd

from pathway_forte.prediction.multiclass import (
    get_class_labels, get_sample_ids_with_cancer_subtypes, match_samples, stabilize_ssgsea_scores_df,
    train_multiclass_svm,
)

__all__ = [
    'do_subtype_prediction',
]

logger = logging.getLogger(__name__)


def do_subtype_prediction(
        ssgsea,
        subtypes,
        *,
        outer_cv_splits,
        inner_cv_splits,
        chain_pca,
        explained_variance,
):
    patient_ids = get_sample_ids_with_cancer_subtypes(subtypes)
    brca_subtypes_df = pd.read_csv(subtypes, sep='\t')

    enrichment_score_df = stabilize_ssgsea_scores_df(ssgsea)
    pathway_features = match_samples(enrichment_score_df, patient_ids)
    class_labels = get_class_labels(pathway_features, brca_subtypes_df)

    all_metrics = train_multiclass_svm(
        pathway_features,
        class_labels,
        outer_cv_splits=outer_cv_splits,
        inner_cv_splits=inner_cv_splits,
        chain_pca=chain_pca,
        explained_variance=explained_variance,
    )
    logger.info(all_metrics)
