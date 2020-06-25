# -*- coding: utf-8 -*-

"""Utilities module."""

import io
import os
import urllib.request
import zipfile
from io import StringIO

import matplotlib.pyplot as plt
import requests
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


def handle_zipfile_download(path: str) -> None:
    """Download and extract zip file content."""
    r = requests.get(path)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall()


def handle_file_download(url: str, filename: str):
    """Handle file download from url, write to file and save to static directory."""
    response = StringIO(urllib.request.urlopen(url).read().decode('utf-8'))
    with open(filename, 'w') as f:
        f.write(response.read())
