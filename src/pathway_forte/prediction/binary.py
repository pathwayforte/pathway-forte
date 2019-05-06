# -*- coding: utf-8 -*-

"""Elastic Net regression with nested cross validation module."""

import gseapy
import logging
import numpy as np
import os
import pandas as pd
from pathway_forte.constants import CLASSIFIER_RESULTS
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from tqdm import tqdm

log = logging.getLogger(__name__)
np.random.seed(10)

__all__ = [
    'ssgsea_nes_to_df',
    'get_parameter_values',
    'train_elastic_net_model',
]


def ssgsea_nes_to_df(ssgsea_scores_csv, classes_file, removed_random=None):
    """Create dataFrame of Normalized Enrichment Scores (NES) from ssGSEA of TCGA expression data.

    :param ssgsea_scores_csv: Text file containing normalized ES for pathways from each sample
    :param test_size: Default test size is 0.25
    :param Optional[int] removed_random: Remove percentage of df
    :return:
    :rtype:
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
    phenoA, phenoB, class_vector = gseapy.parser.gsea_cls_parser(classes_file)

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


def get_parameter_values():
    """Return a list of values that are used by the elastic net as hyperparameters."""
    numbers = [
        i
        for i in range(0, 101)
    ]

    parameters = []

    for i in numbers:

        if (i < 70 and (i % 2) == 0) or ((i % 3) == 0) or ((i % 5) == 0):
            continue

        parameters.append(i / 100)

    return parameters


def train_elastic_net_model(
        x_features,
        y_labels,
        outer_cv_splits,
        inner_cv_splits,
        hyperparameter_space,
        model_name,
        max_iter=1000,
        export=True
):
    """Train elastic net model within a defined hyperparameter space via a nested cross validation given
    expression data.

    :param numpy.array x_features: 2D matrix of pathway scores and samples
    :param list y_labels: class labels of samples
    :param int outer_cv_splits: number of folds for cross validation split in outer loop
    :param int inner_cv_splits: number of folds for cross validation split in inner loop
    :param list hyperparameter_space: list of hyperparameters for l1 and l2 priors
    :param str model_name: name of the model
    :param int max_iter: default to 1000 to ensure convergence
    :return:
    """

    # Use variation of KFold cross validation that returns stratified folds for outer loop in the CV.
    # The folds are made by preserving the percentage of samples for each class.
    skf = StratifiedKFold(n_splits=outer_cv_splits, shuffle=True)

    auc_scores = []

    iterator = tqdm(skf.split(x_features, y_labels))

    for i, (train_index, test_index) in enumerate(iterator):
        X_train, X_test = x_features.iloc[train_index], x_features.iloc[test_index]
        y_train, y_test = [y_labels[i] for i in train_index], [y_labels[i] for i in test_index]

        # Instantiate the model fitting along a regularization path (CV).
        # Inner loop
        glm_elastic = ElasticNetCV(hyperparameter_space, cv=inner_cv_splits, max_iter=max_iter)
        glm_elastic.fit(X_train, y_train)

        log.info('Iteration {}:{}'.format(i, glm_elastic.get_params()))

        predicted_values = glm_elastic.predict(X_test)

        score = roc_auc_score(y_test, predicted_values)
        auc_scores.append(score)

        log.info('Iteration {} AUC score:{}'.format(i, score))

        # Pickle to model
        if export:
            from joblib import dump
            dump(glm_elastic, os.path.join(CLASSIFIER_RESULTS, '{}_{}.joblib'.format(model_name, i)))

    return auc_scores
