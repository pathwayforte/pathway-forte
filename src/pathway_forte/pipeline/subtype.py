# -*- coding: utf-8 -*-

"""CLI wrapper to perform subtype classification."""

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
    # Adapt ssGSEA dataframe for scikit learn purposes
    enrichment_score_df = stabilize_ssgsea_scores_df(ssgsea_path)
    # Prepare sample ids file
    patient_ids = get_sample_ids_with_cancer_subtypes(subtypes_path)
    # Remove indexes from patients that are not in any subtype
    filter_by_index(enrichment_score_df, patient_ids)
    # Read subtype dataframe
    subtypes_df = pd.read_csv(subtypes_path, sep='\t')
    # Prepare class label vector
    class_labels = get_class_labels(enrichment_score_df, subtypes_df)

    # Call the method to train the classifier
    return train_multiclass_classifier(
        enrichment_score_df,
        class_labels,
        outer_cv_splits=outer_cv_splits,
        inner_cv_splits=inner_cv_splits,
        chain_pca=chain_pca,
        explained_variance=explained_variance,
    )
