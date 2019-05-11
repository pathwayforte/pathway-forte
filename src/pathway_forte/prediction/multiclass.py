# -*- coding: utf-8 -*-

"""Logistic regression with nested cross validation module."""

from typing import Optional

import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.multiclass import OneVsOneClassifier
from sklearn.svm import SVC
from tqdm import tqdm

from pathway_forte.prediction.utils import pca_chaining


def get_sample_ids_with_cancer_subtypes(file_path):
    subtype_column = 'subtype_BRCA_Subtype_PAM50'
    df = pd.read_csv(file_path, sep='\t')
    df.head()
    df.drop("barcode", axis=1, inplace=True)  # TODO: remove barcode column upon creation
    df.drop(df[df[subtype_column] == 'Normal'].index, inplace=True)

    return df.index


def stabilize_ssgsea_scores_df(file_path):
    df = pd.read_csv(file_path, sep='\t', header=0)

    # Transpose dataFrame to arrange columns as pathways and rows as genes
    df = df.transpose()

    # Set column index to the first row in the dataframe
    df.columns = df.iloc[0]

    # Remove the first row because it is already set as column index
    df = df.drop("Term|NES")

    return df


def match_samples(df1, list1):
    # Filter out samples with no matches in both datasets
    for index, ssgsea_row in df1.iterrows():
        if index not in list1:
            df1.drop(index, inplace=True)

    return df1


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


PARAM_GRID = {
    'estimator__kernel': ('linear', 'rbf'),
    'estimator__C': [1, 10, 100],
}


def train_multiclass_svm(
        x,
        y,
        inner_cv: int,
        outer_cv: int,
        chain_pca: bool = False,
        explained_variance: Optional[float] = None,
):
    """Train SVM with multiclass labels with a defined hyper-parameter space via a nested cross validation for TCGA
    expression data.

    :param Numpy.ndarray x: 2D array of pathway scores and samples
    :param Numpy.ndarray y: 1D array of sample subtype labels
    :param inner_cv: number of folds for cross validation split in inner loop
    :param outer_cv: number of folds for cross validation split in outer loop
    :param chain_pca: chain PCA to classifier
    :param explained_variance: amount of variance retained. Defaults to 0.95.
    :return:
    """
    all_accuracy_metrics = {}
    all_f1_metrics = {}

    target_names = ['Class 0', 'Class 1', 'Class 2', 'Class 3']

    it = _help_train_multiclass_svm(x, y, inner_cv, outer_cv, chain_pca, explained_variance)
    for i, (classifier, y_test, y_pred) in enumerate(it, start=1):
        # Get the subset accuracy st labels predicted for a sample exactly match true labels (harsh)
        accuracy = metrics.accuracy_score(y_test, y_pred)  # set sample_weight to get weighted accuracy
        f1_score = metrics.f1_score(y_test, y_pred, average="weighted")
        all_accuracy_metrics[i] = accuracy
        all_f1_metrics[i] = f1_score

        print(f'For iteration {i}:')
        print(f'best parameter is {classifier.best_params_}')
        print(f'test accuracy is {accuracy}')
        print(f'f1 score is {f1_score}\n')
        print(metrics.classification_report(y_test, y_pred, target_names=target_names))

    return all_accuracy_metrics, all_f1_metrics


def _help_train_multiclass_svm(
        x,
        y,
        inner_cv: int,
        outer_cv: int,
        chain_pca: bool = False,
        explained_variance: Optional[float] = None,
):
    """Train SVM with multiclass labels with a defined hyper-parameter space via a nested cross validation for TCGA
    expression data.

    :param Numpy.ndarray x: 2D array of pathway scores and samples
    :param Numpy.ndarray y: 1D array of sample subtype labels
    :param inner_cv: number of folds for cross validation split in inner loop
    :param outer_cv: number of folds for cross validation split in outer loop
    :param chain_pca: chain PCA to classifier
    :param explained_variance: amount of variance retained
    :return:
    """
    explained_variance = explained_variance or 0.95

    kf = KFold(n_splits=outer_cv, shuffle=True)

    iterator = tqdm(kf.split(x, y))

    for train_index, test_index in iterator:
        x_train = x[train_index]
        x_test = x[test_index]
        y_train = np.asarray([y[i] for i in train_index])
        y_test = np.asarray([y[i] for i in test_index])

        if chain_pca:  # Apply PCA
            x_train, x_test = pca_chaining(x_train, x_test, explained_variance)

        # Fit one classifier per class
        # For each classifier, class is fit against all other classes
        svm = OneVsOneClassifier(SVC(gamma='scale'))

        classifier = GridSearchCV(estimator=svm, param_grid=PARAM_GRID, cv=inner_cv)
        classifier.fit(x_train, y_train)

        y_pred = classifier.predict(x_test)

        yield classifier, y_test, y_pred
