# -*- coding: utf-8 -*-

"""Logistic regression with nested cross validation module."""

from typing import Collection, Optional

import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.multiclass import OneVsOneClassifier
from sklearn.svm import SVC
from tqdm import tqdm

from pathway_forte.prediction.utils import pca_chaining

SVC_PARAM_GRID = {
    'estimator__kernel': ('linear', 'rbf'),
    'estimator__C': [1, 10, 100],
}


def get_sample_ids_with_cancer_subtypes(path: str, subtype_column: str = 'subtype_BRCA_Subtype_PAM50'):
    df = pd.read_csv(path, sep='\t')
    df.drop("barcode", axis=1, inplace=True)  # TODO: remove barcode column upon creation
    df.drop(df[df[subtype_column] == 'Normal'].index, inplace=True)
    return df.index


def stabilize_ssgsea_scores_df(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', header=0)

    # Transpose dataFrame to arrange columns as pathways and rows as genes
    df = df.transpose()

    # Set column index to the first row in the dataframe
    df.columns = df.iloc[0]

    # Remove the first row because it is already set as column index
    df = df.drop("Term|NES")

    return df


def filter_by_index(df: pd.DataFrame, keep_indexes: Collection) -> None:
    """Filter out samples with no matches in both datasets."""
    for index, ssgsea_row in df.iterrows():
        if index not in keep_indexes:
            df.drop(index, inplace=True)


def get_class_labels(features_df, labels_df):
    # Merge dataFrames and match sample ssGSEA scores with their cancer subtypes
    merged_dfs = pd.merge(features_df, labels_df, left_index=True, right_index=True)

    labels = merged_dfs['subtype_BRCA_Subtype_PAM50'].tolist()

    # change label names to int
    label_dict = {'LumA': 0, 'LumB': 1, 'Her2': 2, 'Basal': 3}
    class_labels = [label_dict[name] for name in labels]
    class_labels = np.asarray(class_labels, dtype=int)

    return class_labels


def convert_df_to_features_array(df):
    # Get list of pathways as features
    feature_cols = list(df.columns.values)

    # Features
    pathways = df[feature_cols]  # Features

    # Transform features dataFrame to numpy array
    pathways_array = pathways.values

    return np.asarray(pathways_array)


def train_multiclass_classifier(
        x,
        y,
        *,
        inner_cv_splits: int,
        outer_cv_splits: int,
        chain_pca: bool = False,
        explained_variance: Optional[float] = None,
        get_estimator=None,
):
    """Train SVM with multiclass labels with a defined hyper-parameter space via a nested cross validation for TCGA
    expression data.

    :param Numpy.ndarray x: 2D array of pathway scores and samples
    :param Numpy.ndarray y: 1D array of sample subtype labels
    :param inner_cv_splits: number of folds for cross validation split in inner loop
    :param outer_cv_splits: number of folds for cross validation split in outer loop
    :param chain_pca: chain PCA to classifier
    :param explained_variance: amount of variance retained. Defaults to 0.95.
    :param get_estimator: A function returning a pair of the classifier and the param grid dictionary to use
     during grid search CV. Defaults to an SVM classifier with :data:`SVC_PARAM_GRID`.
    :return:
    """
    it = _help_train_multiclass_classifier(
        x,
        y,
        inner_cv_splits=inner_cv_splits,
        outer_cv_splits=outer_cv_splits,
        chain_pca=chain_pca,
        explained_variance=explained_variance,
        get_estimator=get_estimator,
    )
    for classifier, y_test, y_pred in it:
        # Get the subset accuracy st labels predicted for a sample exactly match true labels (harsh)
        yield {
            'evaluation': {
                'best_parameters': classifier.best_params_,
                'metrics': {
                    'accuracy': metrics.accuracy_score(y_test, y_pred),  # set sample_weight to get weighted accuracy,
                    'f1': metrics.f1_score(y_test, y_pred, average="weighted"),
                    'mcc': metrics.matthews_corrcoef(y_test, y_pred),
                },
            },
            'data': {
                'y_test': list(y_test),
                'y_pred': list(y_pred),
            },
        }


def get_ovo_svc_classifier():
    """Get a classifier that fits one classifier per class.

    For each classifier, class is fit against all other classes.
    """
    classifier = OneVsOneClassifier(SVC(gamma='scale'))
    return classifier, SVC_PARAM_GRID


def _help_train_multiclass_classifier(
        x,
        y,
        *,
        inner_cv_splits: int,
        outer_cv_splits: int,
        chain_pca: bool = False,
        explained_variance: Optional[float] = None,
        get_estimator=None,
):
    """Train a multi-class classifier over a defined parameter grid using nested cross-validation.

    :param Numpy.ndarray x: 2D array of pathway scores and samples
    :param Numpy.ndarray y: 1D array of sample subtype labels
    :param inner_cv_splits: number of folds for cross validation split in inner loop
    :param outer_cv_splits: number of folds for cross validation split in outer loop
    :param chain_pca: chain PCA to classifier
    :param explained_variance: amount of variance retained
    :param get_estimator: A function returning a pair of the classifier and the param grid dictionary to use
     during grid search CV. Defaults to an SVM classifier with :data:`SVC_PARAM_GRID`.
    :return:
    """
    explained_variance = explained_variance or 0.95
    get_estimator = get_estimator or get_ovo_svc_classifier

    k_fold = KFold(n_splits=outer_cv_splits, shuffle=True)

    iterator = tqdm(k_fold.split(x, y))

    for train_indexes, test_indexes in iterator:
        x_train = x[train_indexes]
        x_test = x[test_indexes]
        y_train = np.asarray([y[train_index] for train_index in train_indexes])
        y_test = np.asarray([y[test_index] for test_index in test_indexes])

        if chain_pca:  # Apply PCA
            x_train, x_test = pca_chaining(x_train, x_test, explained_variance)

        estimator, param_grid = get_estimator()

        classifier = GridSearchCV(estimator=estimator, param_grid=param_grid, cv=inner_cv_splits)
        classifier.fit(x_train, y_train)

        y_pred = classifier.predict(x_test)

        yield classifier, y_test, y_pred
