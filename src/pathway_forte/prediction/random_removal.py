# -*- coding: utf-8 -*-

"""Module that test the stability of prediction methods by removing random pathways."""

from collections import defaultdict

from pathway_forte.prediction.binary import train_elastic_net_model, ssgsea_nes_to_df


def run_stability_iterations(
        percentage_list,
        number_runs,
        ssgsea_scores_path,
        phenotypes_path,
        outer_cv_splits,
        inner_cv_splits,
        hyperparameter_space,
        model_name,
        max_iter=1000,
        export=True
):
    """Run multiple stability tests by removing different percentages of the dataset.

    :param list[int] percentage_list: list with the percentages of pathways to be removed in every iteration
    :param int number_runs: number of runs in each iteration
    :param numpy.array ssgsea_scores_path: 2D matrix of pathway scores and samples
    :param list phenotypes_path: class labels of samples
    :param int outer_cv_splits: number of folds for cross validation split in outer loop
    :param int inner_cv_splits: number of folds for cross validation split in inner loop
    :param list hyperparameter_space: list of hyperparameters for l1 and l2 priors
    :param str model_name: name of the model
    :param int max_iter: default to 1000 to ensure convergence
    :return:
    """
    results = defaultdict(list)

    # Different percentages of the data set that will be removed (e.g., 10% the data set, 20%, etc.)
    for percentage in percentage_list:

        # For each percentage try with multiple runs since every iteration the percentage will be random
        for i in range(0, number_runs):
            remove_n = percentage / ...

            x_features, y_labels = ssgsea_nes_to_df(ssgsea_scores_path, phenotypes_path, remove_random=remove_n)

            i_results = test_stability(
                percentage,
                x_features,
                y_labels,
                outer_cv_splits,
                inner_cv_splits,
                hyperparameter_space,
                model_name,
                max_iter=1000,
                export=True
            )

            results[percentage].append(i_results)

    return results


def test_stability(
        percentage_to_remove,
        x_features,
        y_labels,
        outer_cv_splits,
        inner_cv_splits,
        hyperparameter_space,
        model_name,
        max_iter=1000,
        export=True
):
    """Train elastic net model within a defined hyperparameter space via a nested cross validation after removing a
    percentage of the data set.

    :param int percentage_to_remove: percentage of pathways that will be removed [0,100]
    :param numpy.array x_features: 2D matrix of pathway scores and samples
    :param list y_labels: class labels of samples
    :param int outer_cv_splits: number of folds for cross validation split in outer loop
    :param int inner_cv_splits: number of folds for cross validation split in inner loop
    :param list hyperparameter_space: list of hyperparameters for l1 and l2 priors
    :param str model_name: name of the model
    :param int max_iter: default to 1000 to ensure convergence
    :return:
    """

    # TODO: Remove X% of the x_features

    results = train_elastic_net_model(
        x_features,
        y_labels,
        outer_cv_splits,
        inner_cv_splits,
        hyperparameter_space,
        model_name,
        max_iter=max_iter,
        export=export
    )

    # TODO: Return the features that were removed
