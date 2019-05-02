# -*- coding: utf-8 -*-

"""Utilities module."""

import os

import matplotlib.pyplot as plt
import seaborn as sns


def plot_aucs(data, database, data_set, export_directory, shuffled=False):
    """Render boxplot with AUC values."""
    sns.boxplot(
        y=data,
        showfliers=False,
        palette="Set2"
    ).set_title('{}ROC-AUC scores for {}'.format(
        "Shuffled " if shuffled else "",
        data_set
    ), fontsize=14)

    plt.ylabel('AUC', fontsize=13)
    plt.savefig(os.path.join(export_directory, "{}_{}_aucs.png".format(database, data_set)))


def get_num_samples(samples_file):
    """Return the number of samples.

    :param str samples_file: file path
    :rtype: int
    :return: number of samples
    """
    with open(samples_file, 'r') as num_samples:
        sample_numbers = num_samples.read().replace('\n', '')

    return int(sample_numbers)
