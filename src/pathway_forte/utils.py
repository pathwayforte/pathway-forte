# -*- coding: utf-8 -*-

"""Utilities module."""

import os

import matplotlib.pyplot as plt
import seaborn as sns


def plot_aucs(
        data,
        database_name: str,
        dataset_name: str,
        export_directory: str,
        shuffled: bool = False,
) -> None:
    """Render boxplot with AUC values."""
    sns.boxplot(
        y=data,
        showfliers=False,
        palette="Set2"
    ).set_title('{}ROC-AUC scores for {}'.format(
        "Shuffled " if shuffled else "",
        dataset_name
    ), fontsize=14)

    # Save plot in a given folder
    plt.ylabel('AUC', fontsize=13)
    plt.savefig(os.path.join(export_directory, f"{database_name}_{dataset_name}_aucs.png"))


def get_num_samples(samples_file_path: str) -> int:
    """Return the number of samples.

    :param samples_file_path: file path
    :return: number of samples
    """
    with open(samples_file_path) as file:
        sample_numbers = file.read().replace('\n', '')

    return int(sample_numbers)
