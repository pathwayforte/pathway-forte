# -*- coding: utf-8 -*-

""""""

import logging
from typing import Optional

import pandas as pd

from pathway_forte.prediction.multiclass import (
    filter_by_index, get_class_labels, get_sample_ids_with_cancer_subtypes, stabilize_ssgsea_scores_df,
    train_multiclass_classifier,
)

__all__ = [
    'do_subtype_prediction',
]

logger = logging.getLogger(__name__)


def do_subtype_prediction(
        ssgsea_path: str,
        subtypes_path: str,
        *,
        outer_cv_splits: int,
        inner_cv_splits: int,
        chain_pca: bool,
        explained_variance: Optional[float] = None,
):
    enrichment_score_df = stabilize_ssgsea_scores_df(ssgsea_path)
    patient_ids = get_sample_ids_with_cancer_subtypes(subtypes_path)
    filter_by_index(enrichment_score_df, patient_ids)
    subtypes_df = pd.read_csv(subtypes_path, sep='\t')
    class_labels = get_class_labels(enrichment_score_df, subtypes_df)

    return train_multiclass_classifier(
        enrichment_score_df,
        class_labels,
        outer_cv_splits=outer_cv_splits,
        inner_cv_splits=inner_cv_splits,
        chain_pca=chain_pca,
        explained_variance=explained_variance,
    )
