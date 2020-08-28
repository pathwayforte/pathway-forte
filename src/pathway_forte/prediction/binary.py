# -*- coding: utf-8 -*-

"""Elastic Net regression with nested cross validation module.

This workflow trains an elastic net model for a binary classification task
(e.g., tumor vs. normal patients). The training is conducted using a nested
cross validation approach (the number of cross validation in both loops can
be selected). The model used can be easily changed since most of the models
in `scikit-learn <https://scikit-learn.org/>`_ (the machine learning
library used by this package) required the same input.
"""

import logging
import os
from typing import List, Optional, Tuple, Union

import gseapy
import numpy as np
import pandas as pd
from sklearn import linear_model, model_selection
from sklearn.metrics import roc_auc_score, average_precision_score, auc
from sklearn.model_selection import StratifiedKFold
from tqdm import tqdm

from pathway_forte.constants import CLASSIFIER_RESULTS

__all__ = [
    'ssgsea_nes_to_df',
    'get_l1_ratios',
    'train_elastic_net_model',
]

logger = logging.getLogger(__name__)


def ssgsea_nes_to_df(ssgsea_scores_csv, classes_file, removed_random: Optional[int] = None):
    """Create dataFrame of Normalized Enrichment Scores (NES) from ssGSEA of TCGA expression data.

    :param ssgsea_scores_csv: Text file containing normalized ES for pathways from each sample
    :param test_size: Default test size is 0.25
    :param removed_random: Remove percentage of df
    """
    # Read tsv file using the first row as the header.
    # Important! gseapy.samples.normalized.es.txt contains a descriptive header which here is omitted
    enrichment_score = pd.read_csv(ssgsea_scores_csv, sep='\t', header=0)

    if removed_random:
        previous_df_shape = enrichment_score.shape
        drop_indices = np.random.choice(enrichment_score.index, removed_random, replace=False)
        enrichment_score = enrichment_score.drop(drop_indices)

        # Check size of df
        assert previous_df_shape[0] == enrichment_score.shape[0] + len(drop_indices), 'Problem removing random rows'
        assert previous_df_shape[1] == enrichment_score.shape[1], 'Columns should have not changed after removing rows'

    # Transpose dataFrame to arrange columns as pathways and rows as genes
    enrichment_score_df = enrichment_score.transpose()

    # Set column index to the first row in the dataframe
    enrichment_score_df.columns = enrichment_score_df.iloc[0]

    # Remove the first row because it is already set as column index
    enrichment_score_df = enrichment_score_df.drop("Term|NES")

    # Get class labels
    _, _, class_vector = gseapy.parser.gsea_cls_parser(classes_file)

    class_labels = []

    for label in class_vector:
        if label == 'Normal':
            class_labels.append(0)
        elif label == 'Tumor':
            class_labels.append(1)

    # Get list of pathways as features
    feature_cols = list(enrichment_score_df.columns.values)

    # split dataset into features and target variable (i.e., normal vs tumor sample labels)
    pathways = enrichment_score_df[feature_cols]  # Features
    pathways.reset_index(drop=True, inplace=True)

    # Transform features dataFrame to numpy array
    pathways_array = pathways.to_numpy()

    return pathways_array, class_labels


def get_l1_ratios():
    """Return a list of values that are used by the elastic net as hyperparameters."""
    return [
        i / 100
        for i in range(0, 101)
        if not _skip_index(i)
    ]


def _skip_index(i):
    return (i < 70 and (i % 2) == 0) or ((i % 3) == 0) or ((i % 5) == 0)


def train_elastic_net_model(
        x,
        y,
        outer_cv_splits: int,
        inner_cv_splits: int,
        l1_ratio: List[float],
        model_name: str,
        max_iter: Optional[int] = None,
        export: bool = True,
) -> Tuple[List[float], List[float]]:
    """Train elastic net model via a nested cross validation given expression data.

    Uses a defined hyperparameter space for l1_ratio.

    :param numpy.array x: 2D matrix of pathway scores and samples
    :param list y: class labels of samples
    :param outer_cv_splits: number of folds for cross validation split in outer loop
    :param inner_cv_splits: number of folds for cross validation split in inner loop
    :param l1_ratio: list of hyper-parameters for l1 and l2 priors
    :param model_name: name of the model
    :param max_iter: default to 1000 to ensure convergence
    :param export: Export the models using :mod:`joblib`
    :return: A list of AUC-ROC scores
    """
    auc_scores = []
    auc_pr_scores = []
    it = _help_train_elastic_net_model(
        x=x,
        y=y,
        outer_cv_splits=outer_cv_splits,
        inner_cv_splits=inner_cv_splits,
        l1_ratio=l1_ratio,
        max_iter=max_iter,
    )

    # Iterator to calculate metrics for each CV step
    for i, (glm_elastic, y_test, y_pred) in enumerate(it):
        logger.info(f'Iteration {i}: {glm_elastic.get_params()}')
        auc_scores.append(roc_auc_score(y_test, y_pred))
        auc_pr_scores.append(average_precision_score(y_test, y_pred))

        # Export a pickle the model of the given CV
        if export:
            import joblib
            joblib.dump(glm_elastic, os.path.join(CLASSIFIER_RESULTS, f'{model_name}_{i}.joblib'))

    # Return a list with all AUC/AUC-PR scores for each CV step
    return auc_scores, auc_pr_scores


def _help_train_elastic_net_model(
        x,
        y,
        outer_cv_splits: int,
        inner_cv_splits: int,
        l1_ratio: Union[float, List[float]],
        max_iter: Optional[int] = None,
):
    max_iter = max_iter or 1000

    # Use variation of KFold cross validation that returns stratified folds for outer loop in the CV.
    # The folds are made by preserving the percentage of samples for each class.
    skf = StratifiedKFold(n_splits=outer_cv_splits, shuffle=True)

    # tqdm wrapper to print the current CV state
    iterator = tqdm(skf.split(x, y), desc='Outer CV for binary prediction')

    # Parameter Grid
    param_grid = dict(l1_ratio=l1_ratio)

    # Iterator over the CVs
    for train_indexes, test_indexes in iterator:
        # Splice the entire data set so only the training and test sets for this CV iter are used
        x_train, x_test = x[train_indexes], x[test_indexes]
        y_train = [y[train_index] for train_index in train_indexes]
        y_test = [y[test_index] for test_index in test_indexes]

        # Instantiate the model fitting along a regularization path (CV).
        # Inner loop
        estimator = linear_model.LogisticRegression(penalty='elasticnet', solver='saga', max_iter=max_iter)

        glm_elastic = model_selection.GridSearchCV(
            estimator=estimator,
            param_grid=param_grid,
            cv=inner_cv_splits,
            scoring=roc_auc_score
        )

        # Fit model with train data
        glm_elastic.fit(x_train, y_train)

        # Predict trained model with test data
        y_pred = glm_elastic.predict(x_test)

        # Return model and y test y predict to calculate prediction metrics
        yield glm_elastic, y_test, y_pred
