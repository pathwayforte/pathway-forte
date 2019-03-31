# -*- coding: utf-8 -*-

"""Logistic regression with nested cross validation module."""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def pca_chaining(X_train, X_test, explained_variance):
    """Chain PCA with logistic regression.

    :param pandas.core.frame.DataFrame X_train: Training set to apply dimensionality reduction to
    :param pandas.core.series.Series X_test: Test set to apply dimensionality reduction to
    :param Optional explained_variance: Amount of variance retained
    :return: array-like, shape (n_samples, n_components)
    """
    # Make an instance of the model
    pca = PCA(explained_variance)

    # Fit PCA on the training set only
    pca.fit(X_train)

    # Transform both the training and test set
    X_train = pca.transform(X_train)
    X_test = pca.transform(X_test)

    return X_train, X_test


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
            continue

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
