# -*- coding: utf-8 -*-

"""Command line interface."""

import pickle
import warnings

import click
import pandas as pd
from gseapy.parser import gsea_gmt_parser

from pathway_forte.constants import *
from pathway_forte.export_genesets_to_gmt import create_geneset_df, export_gmt_files, get_all_pathway_genesets
from pathway_forte.mappings import get_equivalent_mappings_dict
from pathway_forte.pathway_enrichment.functional_class import (
    create_cls_file, filter_gene_exp_data, run_gsea, run_ssgsea,
)
from pathway_forte.pathway_enrichment.over_representation import (
    filter_fold_change_fd, perform_hypergeometric_test, read_fold_change_df,
)
from pathway_forte.prediction.binary import get_l1_ratios, ssgsea_nes_to_df, train_elastic_net_model
from pathway_forte.prediction.multiclass import (
    get_class_labels, get_sample_ids_with_cancer_subtypes, match_samples,
    stabilize_ssgsea_scores_df, train_multiclass_svm,
)
from pathway_forte.prediction.survival import run_survival_all_datasets
from pathway_forte.utils import plot_aucs

logger = logging.getLogger(__name__)


@click.group(help='PathwayForte')
def main():
    """Run PathwayForte."""
    logging.basicConfig(level=20, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


@main.command()
def datasets():
    """List the available cancer datasets."""
    click.echo(f"List of data sets used in the paper:")
    for dataset in sorted(CANCER_DATA_SETS):
        click.echo(f'- {dataset}')


"""Export GMT Files"""


@main.command()
def export():
    """Generate gene set files using ComPath."""
    all_pathway_genesets = get_all_pathway_genesets()
    equivalent_mappings_dict = get_equivalent_mappings_dict()
    df = create_geneset_df(all_pathway_genesets, equivalent_mappings_dict)
    export_gmt_files(df)
    click.echo('Done creating files')


"""ORA Analyses"""


@main.group()
def ora():
    """Perform ORA analysis."""


@ora.command()
@click.option('-d', '--genesets', type=click.Path(exists=True), required=False, help='Path to GMT file')
@click.option('-s', '--fold-changes', type=click.Path(exists=True), required=False, help='Path to fold changes file')
@click.option('--no-threshold', is_flag=True, help='Do not apply threshold')
def fisher(genesets, fold_changes, no_threshold):
    """Performs fisher tests enrichment."""

    # Reverse threshold boolean (if "--no-threshold" threshold=False, else threshold=True)
    threshold = not no_threshold

    if threshold:
        click.echo('Filtering out q values > 0.01 according to fdr_bh')

    fc_df = read_fold_change_df(fold_changes)

    significant_genes = filter_fold_change_fd(fc_df)

    click.echo(f'There are a total of {len(significant_genes)} significant genes.')

    # Note that the parser filters out gene sets smaller than 3 and larger than 1000
    gene_sets = gsea_gmt_parser(genesets)

    enriched_pathways = perform_hypergeometric_test(
        significant_genes,
        gene_sets,
        apply_threshold=threshold
    )

    output = os.path.join(os.getcwd(), 'results.pickle')

    # Export dictionary as pickle
    with open(output, 'wb') as file:
        pickle.dump(enriched_pathways, file, protocol=4)

    click.echo(f'Results exported to {output}')


"""FCS Analyses"""


@main.group()
def fcs():
    """List of FCS Analyses."""


@fcs.command()
@click.option('-d', '--data', type=click.Choice(CANCER_DATA_SETS), required=True, help='Name of dataset')
@click.option('-p', '--permutations', type=int, default=100, show_default=True, help='Number of permutations')
def gsea(data, permutations):
    """Run GSEA on TCGA data."""
    click.echo(f'Running GSEA for the {data} dataset')

    make_gsea_export_directories()

    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')
    num_normal_samples = NORMAL_EXPRESSION_SAMPLES.format(data)
    num_tumor_samples = TUMOR_EXPRESSION_SAMPLES.format(data)

    phenotype_path = PHENOTYPE_CLASSES.format(data)

    if not os.path.isfile(phenotype_path):
        logger.info('Creating cls file')
        create_cls_file(gene_exp, num_normal_samples, num_tumor_samples, data)

    kegg_gene_set, reactome_gene_set, wikipathways_gene_set, merge_gene_set = check_gmt_files()
    _ds = [
        ('KEGG', kegg_gene_set, KEGG_GSEA, KEGG_GSEA_TSV),
        ('Reactome', reactome_gene_set, REACTOME_GSEA, REACTOME_GSEA_TSV),
        ('WikiPathways', wikipathways_gene_set, WIKIPATHWAYS_GSEA, WIKIPATHWAYS_GSEA_TSV),
        ('MergeDataset', merge_gene_set, MERGE_GSEA, MERGE_GSEA_TSV),
    ]
    for name, gene_set, output_dir, tsv_fmt in _ds:
        logger.info(f'Running {name}')
        results = run_gsea(
            gene_exp,
            gene_set,
            phenotype_path,
            permutations=permutations,
            output_dir=output_dir,
        )
        results.res2d.to_csv(tsv_fmt.format(data, TODAY), sep='\t')

    logger.info('Done with GSEA analysis')


@fcs.command()
@click.option('-d', '--data', type=click.Choice(CANCER_DATA_SETS), required=True, help='Name of dataset')
def gsea_msig(data):
    """Run GSEA using MSigDB gene sets."""
    click.echo(f'Running GSEA for the {data} dataset')
    make_gsea_export_directories()
    make_ssgsea_export_directories()

    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')
    num_normal_samples = NORMAL_EXPRESSION_SAMPLES.format(data)
    num_tumor_samples = TUMOR_EXPRESSION_SAMPLES.format(data)

    phenotype_path = PHENOTYPE_CLASSES.format(data)

    if not os.path.isfile(phenotype_path):
        logger.info('Creating cls file')
        create_cls_file(gene_exp, num_normal_samples, num_tumor_samples, data)

    logger.info('Running KEGG')
    kegg_gsea_results = run_gsea(
        gene_exp, MSIGDB_KEGG_GENE_SETS, phenotype_path, permutations=100, output_dir=MSIG_GSEA
    )
    kegg_gsea_results.res2d.to_csv(KEGG_MSIG_GSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running Reactome')
    reactome_gsea_results = run_gsea(
        gene_exp, MSIGDB_REACTOME_GENE_SETS, phenotype_path, permutations=100, output_dir=MSIG_GSEA
    )
    reactome_gsea_results.res2d.to_csv(REACTOME_MSIG_GSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running Concatenated Merge')
    conca_merge_gsea_results = run_gsea(
        gene_exp, CONCATENATED_MERGE_GENE_SETS, phenotype_path, permutations=100, output_dir=MERGE_GSEA
    )
    conca_merge_gsea_results.res2d.to_csv(CONCATENATED_MERGE_GSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running KEGG')
    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, MSIGDB_KEGG_GENE_SETS)
    kegg_ssgsea_results = run_ssgsea(filtered_expression_data, MSIGDB_KEGG_GENE_SETS, output_dir=MSIG_SSGSEA)
    kegg_ssgsea_results.res2d.to_csv(KEGG_MSIG_SSGSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running Reactome')
    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, MSIGDB_REACTOME_GENE_SETS)
    reactome_ssgsea_results = run_ssgsea(filtered_expression_data, MSIGDB_REACTOME_GENE_SETS, output_dir=MSIG_SSGSEA)
    reactome_ssgsea_results.res2d.to_csv(REACTOME_MSIG_SSGSEA_TSV.format(data, TODAY), sep='\t')

    logger.info('Running Concatenated Merge')
    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, CONCATENATED_MERGE_GENE_SETS)
    conca_merge_ssgsea_results = run_ssgsea(
        filtered_expression_data, CONCATENATED_MERGE_GENE_SETS, output_dir=MERGE_SSGSEA
    )
    conca_merge_ssgsea_results.res2d.to_csv(CONCATENATED_MERGE_GSEA_TSV.format(data, TODAY), sep='\t')


@fcs.command()
@click.option('-d', '--data', type=click.Choice(CANCER_DATA_SETS), required=True, help='Name of dataset')
def ssgsea(data):
    """Run ssGSEA on TCGA data."""
    click.echo(f'Running ssGSEA for the {data} dataset')

    make_ssgsea_export_directories()

    # Read data
    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')

    kegg_gene_set, reactome_gene_set, wikipathways_gene_set, merge_gene_set = check_gmt_files()

    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, merge_gene_set)

    _ds = [
        ('KEGG', kegg_gene_set, KEGG_SSGSEA, KEGG_SSGSEA_TSV),
        ('Reactome', reactome_gene_set, REACTOME_SSGSEA, REACTOME_SSGSEA_TSV),
        ('WikiPathways', wikipathways_gene_set, WIKIPATHWAYS_SSGSEA, WIKIPATHWAYS_SSGSEA_TSV),
        ('MergeDataset', merge_gene_set, MERGE_SSGSEA, MERGE_SSGSEA_TSV),
    ]
    for name, gene_set, output_dir, fmt in _ds:
        logger.info(f'Running {name}')
        results = run_ssgsea(filtered_expression_data, gene_set, output_dir=output_dir)
        results.res2d.to_csv(fmt.format(data, TODAY), sep='\t')

    logger.info('Done with ssGSEA analysis')


"""Prediction methods."""


@main.group()
def prediction():
    """List of Prediction Methods."""


@prediction.command()
@click.option('-d', '--data', type=click.Choice(CANCER_DATA_SETS), required=True, help='Name of dataset')
@click.option('-ocv', '--outer-cv', type=int, default=10, show_default=True, help='Number of k splits in outer cv')
@click.option('-icv', '--inner-cv', type=int, default=10, show_default=True, help='Number of k splits in inner cv')
@click.option(
    '-i', '--max_iterations', type=int, default=1000, show_default=True,
    help='Number of max iterations to converge'
)
@click.option('--turn-off-warnings', is_flag=True, help='Turns off warnings')
def binary(data, outer_cv, inner_cv, max_iterations, turn_off_warnings):
    """Train elastic net for binary prediction."""
    make_classifier_results_directory()
    click.echo(f'Training Elastic Net via nested CV in the {data} dataset with'
               f' {outer_cv} outer loops and {inner_cv} inner loops')

    if turn_off_warnings:
        click.echo("Warnings are now turned off")
        warnings.simplefilter('ignore')

    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')
    num_normal_samples = NORMAL_EXPRESSION_SAMPLES.format(data)
    num_tumor_samples = TUMOR_EXPRESSION_SAMPLES.format(data)

    if not os.path.isfile(PHENOTYPE_CLASSES.format(data)):
        click.echo('Creating cls file')
        create_cls_file(gene_exp, num_normal_samples, num_tumor_samples, data)

    phenotypes = PHENOTYPE_CLASSES.format(data)

    l1_ratio = get_l1_ratios()
    click.echo(f'L1 ratios: {l1_ratio}')

    kegg_ssgsea_nes_path = os.path.join(KEGG_SSGSEA, f'kegg_{data}.tsv')
    reactome_ssgsea_nes_path = os.path.join(REACTOME_SSGSEA, f'reactome_{data}.tsv')
    wikipathways_ssgsea_nes_path = os.path.join(WIKIPATHWAYS_SSGSEA, f'wikipathways_{data}.tsv')
    merge_ssgsea_nes_path = os.path.join(MERGE_SSGSEA, f'merge_{data}.tsv')
    _ds = [
        ('KEGG', kegg_ssgsea_nes_path, KEGG),
        ('Reactome', reactome_ssgsea_nes_path, REACTOME),
        ('WikiPathways', wikipathways_ssgsea_nes_path, WIKIPATHWAYS),
        ('Merge', merge_ssgsea_nes_path, 'MERGE')
    ]
    for name, tsv_path, fmt in _ds:
        click.echo(f'Training on {name}')
        x, y = ssgsea_nes_to_df(tsv_path, phenotypes)
        aucs = train_elastic_net_model(
            x,
            y,
            outer_cv,
            inner_cv,
            l1_ratio,
            f"{data}-{fmt}",
            max_iterations,
        )
        plot_aucs(aucs, data, fmt, CLASSIFIER_RESULTS)
        click.echo(f"{name} AUCs: {aucs}")


@prediction.command()
@click.option('-d', '--data', type=str, required=True, help='Name of dataset')
@click.option('-ocv', '--outer-cv', type=int, default=10, show_default=True, help='Number of k splits in outer cv')
@click.option('-icv', '--inner-cv', type=int, default=10, show_default=True, help='Number of k splits in inner cv')
@click.option('--turn-off-warnings', is_flag=True, help='Turns off warnings')
def survival(data, outer_cv, inner_cv, turn_off_warnings):
    """Train survival model."""
    click.echo(f'Running survival analysis for {data} with: {outer_cv} outer CVs and {inner_cv} inner CVS')

    param_grid = {
        'l1_ratio': get_l1_ratios(),
    }
    click.echo(f'Parameter grid: {param_grid}')

    if turn_off_warnings:
        click.echo("Warnings are now turned off")
        warnings.simplefilter('ignore')

    results = run_survival_all_datasets(os.path.join(SSGSEA, data), outer_cv, inner_cv, param_grid)

    click.echo(results)

    results_path = os.path.join(CLASSIFIER_RESULTS, f'survival_results_{data}.pickle')

    with open(results_path, 'wb') as file:
        pickle.dump(results, file)

    click.echo(f'Done with survival analysis. Results are exported in {results_path}')


@prediction.command()
@click.option('-d', '--ssgsea', type=click.Path(exists=True), required=True, help='Path to ssGSEA file')
@click.option('-s', '--subtypes', type=click.Path(exists=True), required=True, help='Path to the subtypes file')
@click.option('-ocv', '--outer-cv', type=int, default=10, show_default=True, help='Number of k splits in outer cv')
@click.option('-icv', '--inner-cv', type=int, default=10, show_default=True, help='Number of k splits in inner cv')
@click.option('--turn-off-warnings', is_flag=True, help='Turns off warnings')
def subtype(ssgsea, subtypes, outer_cv, inner_cv, turn_off_warnings):
    """Train subtype analysis."""
    click.echo(f'Running subtype analysis for {ssgsea} with: {outer_cv} outer CVs and {inner_cv} inner CVS')

    if turn_off_warnings:
        click.echo("Warnings are now turned off")
        warnings.simplefilter('ignore')

    patient_ids = get_sample_ids_with_cancer_subtypes(subtypes)
    brca_subtypes_df = pd.read_csv(subtypes, sep='\t')

    enrichment_score_df = stabilize_ssgsea_scores_df(ssgsea)
    pathway_features = match_samples(enrichment_score_df, patient_ids)
    class_labels = get_class_labels(pathway_features, brca_subtypes_df)

    all_metrics = train_multiclass_svm(
        pathway_features,
        class_labels,
        inner_cv=5,
        outer_cv=5,
        chain_pca=False,
        explained_variance=0.95,
    )
    click.echo(all_metrics)


@prediction.command()
@click.option('-s', '--ssgsea-scores-path', type=click.Path(exists=True), required=True,
              help='ssGSEA scores file')
@click.option('-p', '--phenotypes-path', type=click.Path(exists=True), required=True,
              help='Path to the phenotypes file')
@click.option('-ocv', '--outer-cv', type=int, default=10, show_default=True, help='Number of k splits in outer cv')
@click.option('-icv', '--inner-cv', type=int, default=10, show_default=True, help='Number of k splits in inner cv')
@click.option(
    '-i', '--max_iterations', type=int, default=1000, show_default=True,
    help='Number of max iterations to converge'
)
@click.option('--turn-off-warnings', is_flag=True, help='Turns off warnings')
def test_stability_prediction(
        ssgsea_scores_path,
        phenotypes_path,
        outer_cv,
        inner_cv,
        max_iterations,
        turn_off_warnings,
):
    """Train elastic net."""
    make_classifier_results_directory()
    click.echo(
        f'Training Elastic Net via nested CV for {outer_cv} dataset '
        f'with {ssgsea_scores_path} (phenotypes: {phenotypes_path}) '
        f'outer loops and {inner_cv} inner loops'
    )

    if turn_off_warnings:
        click.echo("ssgsea_nes_to_dfWarnings are now turned off")
        warnings.simplefilter('ignore')

    l1_ratio = get_l1_ratios()
    click.echo(f'L1 ratios: {l1_ratio}')

    results = train_elastic_net_model(
        ssgsea_scores_path,
        phenotypes_path,
        outer_cv,
        inner_cv,
        l1_ratio,
        'elastic_net',
        max_iter=max_iterations,
    )

    click.echo(results)


if __name__ == '__main__':
    main()
