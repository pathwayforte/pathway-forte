# -*- coding: utf-8 -*-

"""Utilities for prediction."""

from typing import Tuple

from sklearn.decomposition import PCA

__all__ = [
    'pca_chaining',
]


def pca_chaining(train, test, n_components) -> Tuple:
    """Chain PCA with logistic regression.

    :param pandas.core.frame.DataFrame train: Training set to apply dimensionality reduction to
    :param pandas.core.series.Series test: Test set to apply dimensionality reduction to
    :param n_components: Amount of variance retained
    :return: array-like, shape (n_samples, n_components)
    """
    # Make an instance of the model
    pca = PCA(n_components)

    # Fit PCA on the training set only then transform both
    train_transformed = pca.fit_transform(train)
    test_transformed = pca.transform(test)

    return train_transformed, test_transformed
