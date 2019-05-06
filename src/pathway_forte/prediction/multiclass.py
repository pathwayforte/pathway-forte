# -*- coding: utf-8 -*-

"""Logistic regression with nested cross validation module."""

import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.multiclass import OneVsOneClassifier
from sklearn.svm import SVC
from tqdm import tqdm


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


def convert_df_to_features_array(df):
    # Get list of pathways as features
    feature_cols = list(df.columns.values)

    # Features
    pathways = df[feature_cols]  # Features

    # Transform features dataFrame to numpy array
    pathways_array = pathways.values

    return np.asarray(pathways_array)


def train_multiclass_svm(X, y, inner_cv, outer_cv, chain_pca=False, explained_variance=0.95):
    """Train SVM with multiclass labels with a defined hyperparameter space via a nested cross validation for TCGA
    expression data.

    :param Numpy.ndarray X: 2D array of pathway scores and samples
    :param Numpy.ndarray y: 1D array of sample subtype labels
    :param int inner_cv: number of folds for cross validation split in inner loop
    :param int outer_cv: number of folds for cross validation split in outer loop
    :param chain_pca: chain PCA to classifier
    :param explained_variance: amount of variance retained
    :return:
    """
    all_accuracy_metrics = {}
    all_f1_metrics = {}

    target_names = ['Class 0', 'Class 1', 'Class 2', 'Class 3']

    kf = KFold(n_splits=outer_cv, shuffle=True)

    iterator = tqdm(kf.split(X, y))

    for i, (train_index, test_index) in enumerate(iterator):

        X_train = X[train_index]
        X_test = X[test_index]
        y_train, y_test = np.asarray([y[i] for i in train_index]), np.asarray(
            [y[i] for i in test_index])

        if chain_pca:
            # Apply PCA
            X_train, X_test = pca_chaining(X_train, X_test, explained_variance)

        # Fit one classifier per class
        # For each classifier, class is fit against all other classes
        svm = OneVsOneClassifier(SVC(gamma='scale'))

        # Set up possible values of parameters to optimize over
        p_grid = {'estimator__kernel': ('linear', 'rbf'), 'estimator__C': [1, 10, 100]}

        classifier = GridSearchCV(estimator=svm, param_grid=p_grid, cv=inner_cv)

        classifier.fit(X_train, y_train)
        y_pred = classifier.predict(X_test)

        # Get the subset accuracy st labels predicted for a sample exactly match true labels (harsh)
        accurcay = metrics.accuracy_score(y_test, y_pred)  # set sample_weight to get weighted accuracy
        f1_score = metrics.f1_score(y_test, y_pred, average="weighted")
        all_accuracy_metrics[i + 1] = accurcay
        all_f1_metrics[i + 1] = f1_score

        print('For iteration {}:'.format(i + 1))
        print('best parameter is {}'.format(classifier.best_params_))
        print('test accuracy is {}'.format(accurcay))
        print('f1 score is {}'.format(f1_score))
        print("\n")
        print(metrics.classification_report(y_test, y_pred, target_names=target_names))

    return all_accuracy_metrics, all_f1_metrics
