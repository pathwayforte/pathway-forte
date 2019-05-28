# -*- coding: utf-8 -*-

"""Command line interface."""

import json
import logging
import warnings

import click
from sklearn import metrics

from pathway_forte.constants import CANCER_DATA_SETS
from pathway_forte.pipeline import (
    do_binary_prediction, do_export, do_gsea, do_gsea_msig, do_hypergeometric, do_ssgsea, do_stability_prediction,
    do_subtype_prediction, do_survival_prediction,
)

logger = logging.getLogger(__name__)


@click.group()
def main():
    """Run PathwayForte."""
    logging.basicConfig(level=20, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


@main.command()
def datasets():
    """List the available cancer datasets."""
    click.echo(f"List of data sets used in the paper:")
    for dataset in sorted(CANCER_DATA_SETS):
        click.echo(f'- {dataset}')


@main.command()
def export():
    """Generate gene set files using ComPath."""
    do_export()


@main.group()
def ora():
    """Perform ORA analysis."""


@ora.command()
@click.option('-d', '--genesets', type=click.Path(file_okay=True, dir_okay=False, exists=True), required=True,
              help='Path to GMT file')
@click.option('-s', '--fold-changes', type=click.Path(file_okay=True, dir_okay=False, exists=True), required=True,
              help='Path to fold changes file')
@click.option('--no-threshold', is_flag=True, help='Do not apply threshold')
@click.option('-o', '--output', type=click.Path(), help='Optional path for output JSON file')
def hypergeometric(genesets, fold_changes, no_threshold, output):
    """Performs one-tailed hypergeometric test enrichment."""
    # Reverse threshold boolean (if "--no-threshold" threshold=False, else threshold=True)
    threshold = not no_threshold
    do_hypergeometric(genesets, fold_changes, threshold, output)


@main.group()
def fcs():
    """List of FCS Analyses."""


cancer_data_set_option = click.option(
    '-d', '--data',
    type=click.Choice(CANCER_DATA_SETS),
    required=True,
    help='Name of the cancer dataset from TCGA',
)

outer_cv_option = click.option(
    '--outer-cv',
    type=int,
    default=10,
    show_default=True,
    help='Number of splits in outer cross-validation',
)

inner_cv_option = click.option(
    '--inner-cv',
    type=int,
    default=10,
    show_default=True,
    help='Number of splits in inner cross-validation',
)

turn_off_warnings_option = click.option('--turn-off-warnings', is_flag=True, help='Turns off warnings')


@fcs.command()
@cancer_data_set_option
@click.option('-p', '--permutations', type=int, default=100, show_default=True, help='Number of permutations')
def gsea(data, permutations):
    """Run GSEA on TCGA data."""
    click.echo(f'Running GSEA for the {data} dataset')
    do_gsea(data, permutations)
    click.echo('Done with GSEA analysis')


@fcs.command()
@cancer_data_set_option
def gsea_msig(data):
    """Run GSEA on TCGA data using MSigDB gene sets."""
    click.echo(f'Running GSEA (MSigDB) for the {data} dataset')
    do_gsea_msig(data)
    click.echo(f'Done with GSEA (MSigDB) analysis')


@fcs.command()
@cancer_data_set_option
def ssgsea(data):
    """Run ssGSEA on TCGA data."""
    click.echo(f'Running ssGSEA for the {data} dataset')
    do_ssgsea(data)
    click.echo('Done with ssGSEA analysis')


@main.group()
def prediction():
    """List of Prediction Methods."""


@prediction.command()
@cancer_data_set_option
@outer_cv_option
@inner_cv_option
@click.option(
    '-i', '--max_iterations', type=int, default=1000, show_default=True,
    help='Number of max iterations to converge'
)
@turn_off_warnings_option
def binary(data, outer_cv, inner_cv, max_iterations, turn_off_warnings):
    """Train elastic net for binary prediction."""
    if turn_off_warnings:
        click.echo("Warnings are now turned off")
        warnings.simplefilter('ignore')

    do_binary_prediction(
        data=data,
        outer_cv_splits=outer_cv,
        inner_cv_splits=inner_cv,
        max_iter=max_iterations,
    )


@prediction.command()
@click.option('-d', '--data', type=str, required=True, help='Name of dataset')
@outer_cv_option
@inner_cv_option
@turn_off_warnings_option
def survival(data, outer_cv, inner_cv, turn_off_warnings):
    """Train survival model."""
    if turn_off_warnings:
        click.echo("Warnings are now turned off")
        warnings.simplefilter('ignore')

    click.echo(f'Running survival analysis for {data} with: {outer_cv} outer CVs and {inner_cv} inner CVs')
    do_survival_prediction(data=data, outer_cv_splits=outer_cv, inner_cv_splits=inner_cv)


@prediction.command()
@click.option('-d', '--ssgsea', type=click.Path(exists=True), required=True, help='Path to ssGSEA file')
@click.option('-s', '--subtypes', type=click.Path(exists=True), required=True, help='Path to the subtypes file')
@outer_cv_option
@inner_cv_option
@click.option('--chain-pca', is_flag=True)
@click.option('--explained-variance', type=float, default=0.95, show_default=True, help='Explained variance')
@turn_off_warnings_option
def subtype(
        ssgsea,
        subtypes,
        outer_cv: int,
        inner_cv: int,
        explained_variance: float,
        chain_pca: bool,
        turn_off_warnings: bool,
):
    """Train subtype analysis."""
    if turn_off_warnings:
        click.echo("Warnings are now turned off")
        warnings.simplefilter('ignore')

    click.echo(f'Running subtype analysis for {ssgsea} with: {outer_cv} outer CVs and {inner_cv} inner CVs')
    results = do_subtype_prediction(
        ssgsea,
        subtypes,
        outer_cv_splits=outer_cv,
        inner_cv_splits=inner_cv,
        chain_pca=chain_pca,
        explained_variance=explained_variance,
    )
    for i, result in enumerate(results, start=1):
        click.echo(json.dumps(result['evaluation'], indent=2))
        click.echo(metrics.classification_report(
            result['data']['y_test'],
            result['data']['y_pred'],
            target_names=['Class 0', 'Class 1', 'Class 2', 'Class 3'],
        ))


@prediction.command()
@click.option('-s', '--ssgsea-scores-path', type=click.Path(exists=True), required=True,
              help='ssGSEA scores file')
@click.option('-p', '--phenotypes-path', type=click.Path(exists=True), required=True,
              help='Path to the phenotypes file')
@outer_cv_option
@inner_cv_option
@click.option(
    '-i', '--max_iterations', type=int, default=1000, show_default=True,
    help='Number of max iterations to converge'
)
@turn_off_warnings_option
def test_stability_prediction(
        ssgsea_scores_path,
        phenotypes_path,
        outer_cv,
        inner_cv,
        max_iterations,
        turn_off_warnings,
):
    if turn_off_warnings:
        click.echo("ssgsea_nes_to_dfWarnings are now turned off")
        warnings.simplefilter('ignore')

    do_stability_prediction(
        ssgsea_scores_path=ssgsea_scores_path,
        phenotypes_path=phenotypes_path,
        outer_cv_splits=outer_cv,
        inner_cv_splits=inner_cv,
        max_iter=max_iterations,
    )


if __name__ == '__main__':
    main()
