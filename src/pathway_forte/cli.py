# -*- coding: utf-8 -*-

"""Command line interface."""

import pickle
import warnings

import click
from gseapy.parser import gsea_gmt_parser

from pathway_forte.export_genesets_to_gmt import *
from pathway_forte.mappings import *
from pathway_forte.pathway_enrichment.functional_class import (
    create_cls_file,
    run_gsea,
    run_ssgsea,
    filter_gene_exp_data,
)
from pathway_forte.pathway_enrichment.over_representation import (
    read_fold_change_df,
    filter_fold_change_fd,
    perform_hypergeometric_test
)
from pathway_forte.prediction.binary import get_parameter_values, ssgsea_nes_to_df, train_elastic_net_model
from pathway_forte.prediction.multiclass import *
from pathway_forte.prediction.survival import run_survival_all_datasets
from pathway_forte.utils import plot_aucs

logger = logging.getLogger(__name__)

date_today = time.strftime("%d_%m_%Y")


@click.group(help='PathwayForte')
def main():
    """Run PathwayForte."""
    logging.basicConfig(level=20, format="%(asctime)s - %(levelname)s - %(name)s - %(message)s")


@main.command()
def datasets():
    """List of Cancer Datasets."""
    click.echo("List of data sets used in the paper: {}".format(CANCER_DATA_SETS))


"""Export GMT Files"""


@main.command()
def export_gene_sets():
    """Generate Gene Sets using ComPath."""
    from bio2bel_kegg import Manager as KeggManager
    from bio2bel_reactome import Manager as ReactomeManager
    from bio2bel_wikipathways import Manager as WikipathwaysManager
    """Export GMT files."""
    # Initialize managers
    kegg_manager = KeggManager()
    # Initiate Reactome Manager
    reactome_manager = ReactomeManager()
    # Initiate WikiPathways Manager
    wikipathways_manager = WikipathwaysManager()

    kegg_reactome_df, kegg_wikipathways_df, wikipathways_reactome_df, special_mappings_df = load_compath_mapping_dfs()

    # Get mapping dictionary of all equivalent pairs of pathways. Special mappings are not included in the overall
    # mappings as some of the WP pathways possess identical IDs.

    equivalent_mappings_dict = get_mapping_dict(
        pd.concat([kegg_reactome_df, kegg_wikipathways_df, wikipathways_reactome_df]),
        'equivalentTo'
    )

    logger.info('Getting ComPath/Bio2BEL KEGG gene sets')
    kegg_gene_set = get_compath_genesets(KEGG, kegg_manager)

    logger.info('Getting ComPath/Bio2BEL Reactome gene sets')
    reactome_gene_set = get_compath_genesets(REACTOME, reactome_manager)

    logger.info('Getting ComPath/Bio2BEL WikiPathways gene sets')
    wikipathways_gene_set = get_compath_genesets(WIKIPATHWAYS, wikipathways_manager)

    # Get ComPath/Bio2BEL genesets for each pathway, from each database
    all_pathway_genesets = {**kegg_gene_set, **reactome_gene_set, **wikipathways_gene_set}

    df = create_geneset_df(all_pathway_genesets, equivalent_mappings_dict)

    export_gmt_files(df)
    click.echo('Done creating files')


"""ORA Analyses"""


@main.group()
def ora():
    """List of ORA Analyses."""


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
    click.echo('Running GSEA for the {} dataset'.format(data))

    make_gsea_export_directories()

    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')
    num_normal_samples = NORMAL_EXPRESSION_SAMPLES.format(data)
    num_tumor_samples = TUMOR_EXPRESSION_SAMPLES.format(data)

    kegg_gene_set, reactome_gene_set, wikipathways_gene_set, merge_gene_set = check_gmt_files()

    if not os.path.isfile(PHENOTYPE_CLASSES.format(data)):
        logger.info('Creating cls file')
        create_cls_file(gene_exp, num_normal_samples, num_tumor_samples, data)

    PHENOTYPE = PHENOTYPE_CLASSES.format(data)

    logger.info('Running KEGG')
    kegg_gsea_results = run_gsea(
        gene_exp, kegg_gene_set, PHENOTYPE, permutations=permutations, output_dir=KEGG_GSEA
    )

    kegg_gsea_results.res2d.to_csv(KEGG_GSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running Reactome')
    reactome_gsea_results = run_gsea(
        gene_exp, reactome_gene_set, PHENOTYPE, permutations=permutations, output_dir=REACTOME_GSEA
    )

    reactome_gsea_results.res2d.to_csv(REACTOME_GSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running WikiPathways')
    wikipathways_gsea_results = run_gsea(
        gene_exp, wikipathways_gene_set, PHENOTYPE, permutations=permutations, output_dir=WIKIPATHWAYS_GSEA
    )

    wikipathways_gsea_results.res2d.to_csv(WIKIPATHWAYS_GSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running MergeDataset')
    merge_gsea_results = run_gsea(
        gene_exp, merge_gene_set, PHENOTYPE, permutations=permutations, output_dir=MERGE_GSEA
    )

    merge_gsea_results.res2d.to_csv(MERGE_GSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Done with GSEA analysis')


@fcs.command()
@click.option('-d', '--data', type=click.Choice(CANCER_DATA_SETS), required=True, help='Name of dataset')
def gsea_msig(data):
    """Run GSEA using MSigDB gene sets."""
    click.echo('Running GSEA for the {} dataset'.format(data))
    make_gsea_export_directories()
    make_ssgsea_export_directories()

    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')
    num_normal_samples = NORMAL_EXPRESSION_SAMPLES.format(data)
    num_tumor_samples = TUMOR_EXPRESSION_SAMPLES.format(data)

    if not os.path.isfile(PHENOTYPE_CLASSES.format(data)):
        logger.info('Creating cls file')
        create_cls_file(gene_exp, num_normal_samples, num_tumor_samples, data)

    PHENOTYPE = PHENOTYPE_CLASSES.format(data)

    logger.info('Running KEGG')
    kegg_gsea_results = run_gsea(
        gene_exp, MSIGDB_KEGG_GENE_SETS, PHENOTYPE, permutations=100, output_dir=MSIG_GSEA
    )

    kegg_gsea_results.res2d.to_csv(KEGG_MSIG_GSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running Reactome')
    reactome_gsea_results = run_gsea(
        gene_exp, MSIGDB_REACTOME_GENE_SETS, PHENOTYPE, permutations=100, output_dir=MSIG_GSEA
    )
    reactome_gsea_results.res2d.to_csv(REACTOME_MSIG_GSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running Concatenated Merge')
    conca_merge_gsea_results = run_gsea(
        gene_exp, CONCATENATED_MERGE_GENE_SETS, PHENOTYPE, permutations=100, output_dir=MERGE_GSEA
    )

    conca_merge_gsea_results.res2d.to_csv(CONCATENATED_MERGE_GSEA_TSV.format(data, date_today), sep='\t')

    click.echo('Running ssGSEA for the {} dataset'.format(data))

    logger.info('Running KEGG')

    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, MSIGDB_KEGG_GENE_SETS)
    kegg_ssgsea_results = run_ssgsea(filtered_expression_data, MSIGDB_KEGG_GENE_SETS, output_dir=MSIG_SSGSEA)

    kegg_ssgsea_results.res2d.to_csv(KEGG_MSIG_SSGSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running Reactome')

    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, MSIGDB_REACTOME_GENE_SETS)
    reactome_ssgsea_results = run_ssgsea(filtered_expression_data, MSIGDB_REACTOME_GENE_SETS, output_dir=MSIG_SSGSEA)

    reactome_ssgsea_results.res2d.to_csv(REACTOME_MSIG_SSGSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running Concatenated Merge')

    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, CONCATENATED_MERGE_GENE_SETS)
    conca_merge_ssgsea_results = run_ssgsea(
        filtered_expression_data, CONCATENATED_MERGE_GENE_SETS, output_dir=MERGE_SSGSEA
    )

    conca_merge_ssgsea_results.res2d.to_csv(CONCATENATED_MERGE_GSEA_TSV.format(data, date_today), sep='\t')


@fcs.command()
@click.option('-d', '--data', type=click.Choice(CANCER_DATA_SETS), required=True, help='Name of dataset')
def ssgsea(data):
    """Run ssGSEA on TCGA data."""
    click.echo('Running ssGSEA for the {} dataset'.format(data))

    make_ssgsea_export_directories()

    kegg_gene_set, reactome_gene_set, wikipathways_gene_set, merge_gene_set = check_gmt_files()

    # Read data
    gene_exp = pd.read_csv(EXPRESSION_MATRIX.format(data), sep='\t')

    # Filter data set such that only those genes which are in the gene sets are in the expression data
    filtered_expression_data = filter_gene_exp_data(gene_exp, merge_gene_set)

    logger.info('Running KEGG')
    kegg_ssgsea_results = run_ssgsea(filtered_expression_data, kegg_gene_set, output_dir=KEGG_SSGSEA)

    kegg_ssgsea_results.res2d.to_csv(KEGG_SSGSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running Reactome')
    reactome_ssgsea_results = run_ssgsea(filtered_expression_data, reactome_gene_set, output_dir=REACTOME_SSGSEA)

    reactome_ssgsea_results.res2d.to_csv(REACTOME_SSGSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running WikiPathways')
    wikipathways_ssgsea_results = run_ssgsea(
        filtered_expression_data, wikipathways_gene_set, output_dir=WIKIPATHWAYS_SSGSEA
    )

    wikipathways_ssgsea_results.res2d.to_csv(WIKIPATHWAYS_SSGSEA_TSV.format(data, date_today), sep='\t')

    logger.info('Running MergeDataset')
    merge_ssgsea_results = run_ssgsea(filtered_expression_data, merge_gene_set, output_dir=MERGE_SSGSEA)

    merge_ssgsea_results.res2d.to_csv(MERGE_SSGSEA_TSV.format(data, date_today), sep='\t')
    #
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
    click.echo('Training Elastic Net via nested CV in the {} dataset with {} outer loops and {} inner loops'.format(
        data,
        outer_cv,
        inner_cv
    ))

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

    KEGG_SSGSEA_NES = os.path.join(KEGG_SSGSEA, 'kegg_{}.tsv'.format(data))
    REACTOME_SSGSEA_NES = os.path.join(REACTOME_SSGSEA, 'reactome_{}.tsv'.format(data))
    WIKIPATHWAYS_SSGSEA_NES = os.path.join(WIKIPATHWAYS_SSGSEA, 'wikipathways_{}.tsv'.format(data))
    MERGE_SSGSEA_NES = os.path.join(MERGE_SSGSEA, 'merge_{}.tsv'.format(data))

    parameter_list = get_parameter_values()
    click.echo('Hyperparameter list {}'.format(parameter_list))

    """KEGG"""
    click.echo('Training on KEGG')

    x, y = ssgsea_nes_to_df(KEGG_SSGSEA_NES, phenotypes)

    aucs = train_elastic_net_model(x, y, outer_cv, inner_cv, parameter_list, "{}-{}".format(data, KEGG), max_iterations)

    plot_aucs(aucs, data, KEGG, CLASSIFIER_RESULTS)

    click.echo("KEGG AUCS: {}".format(aucs))

    """Reactome"""
    click.echo('Training on Reactome')

    x, y = ssgsea_nes_to_df(REACTOME_SSGSEA_NES, phenotypes)

    aucs = train_elastic_net_model(
        x, y, outer_cv, inner_cv, parameter_list, "{}-{}".format(data, REACTOME), max_iterations
    )

    click.echo("Reactome AUCS: {}".format(aucs))

    plot_aucs(aucs, data, REACTOME, CLASSIFIER_RESULTS)

    """WikiPathways"""
    click.echo('Training on WikiPathways')

    x, y = ssgsea_nes_to_df(WIKIPATHWAYS_SSGSEA_NES, phenotypes)

    aucs = train_elastic_net_model(
        x, y, outer_cv, inner_cv, parameter_list, "{}-{}".format(data, WIKIPATHWAYS), max_iterations
    )

    click.echo("WikiPathways AUCS: {}".format(aucs))

    plot_aucs(aucs, data, WIKIPATHWAYS, CLASSIFIER_RESULTS)

    """Merge"""
    click.echo('Training on Merge')
    x, y = ssgsea_nes_to_df(MERGE_SSGSEA_NES, phenotypes)

    aucs = train_elastic_net_model(
        x, y, outer_cv, inner_cv, parameter_list, "{}-{}".format(data, "MERGE"), max_iterations
    )

    plot_aucs(aucs, data, 'MERGE', CLASSIFIER_RESULTS)

    click.echo("MERGE_SSGSEA_NES AUCS: {}".format(aucs))


@prediction.command()
@click.option('-d', '--data', type=str, required=True, help='Name of dataset')
@click.option('-ocv', '--outer-cv', type=int, default=10, show_default=True, help='Number of k splits in outer cv')
@click.option('-icv', '--inner-cv', type=int, default=10, show_default=True, help='Number of k splits in inner cv')
@click.option('--turn-off-warnings', is_flag=True, help='Turns off warnings')
def survival(data, outer_cv, inner_cv, turn_off_warnings):
    """Train survival model."""
    click.echo('Running survival analysis for {} with: {} outer CVs and {} inner CVS'.format(data, outer_cv, inner_cv))

    parameter_list = {'l1_ratio': get_parameter_values()}
    click.echo('Hyperparameter list {}'.format(parameter_list))

    if turn_off_warnings:
        click.echo("Warnings are now turned off")
        warnings.simplefilter('ignore')

    results = run_survival_all_datasets(os.path.join(SSGSEA, data), outer_cv, inner_cv, parameter_list)

    click.echo(results)

    results_path = os.path.join(CLASSIFIER_RESULTS, 'survival_results_{}.pickle'.format(data))

    with open(results_path, 'wb') as file:
        pickle.dump(results, file)

    click.echo('Done with survival analysis. Results are exported in {}'.format(results_path))


@prediction.command()
@click.option('-d', '--ssgsea', type=click.Path(exists=True), required=True, help='Path to ssGSEA file')
@click.option('-s', '--subtypes', type=click.Path(exists=True), required=True, help='Path to the subtypes file')
@click.option('-ocv', '--outer-cv', type=int, default=10, show_default=True, help='Number of k splits in outer cv')
@click.option('-icv', '--inner-cv', type=int, default=10, show_default=True, help='Number of k splits in inner cv')
@click.option('--turn-off-warnings', is_flag=True, help='Turns off warnings')
def subtype(ssgsea, subtypes, outer_cv, inner_cv, turn_off_warnings):
    """Train subtype analysis."""
    click.echo('Running subtype analysis for {} with: {} outer CVs and {} inner CVS'.format(ssgsea, outer_cv, inner_cv))

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
        explained_variance=0.95
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
def test_stability_prediction(ssgsea_scores_path, phenotypes_path, outer_cv, inner_cv, max_iterations,
                              turn_off_warnings):
    """Train elastic net."""
    make_classifier_results_directory()
    click.echo(
        'Training Elastic Net via nested CV for {} dataset with {} (phenotypes: {}) outer loops and {} inner loops'.format(
            outer_cv,
            ssgsea_scores_path,
            phenotypes_path,
            inner_cv,

        ))

    if turn_off_warnings:
        click.echo("ssgsea_nes_to_dfWarnings are now turned off")
        warnings.simplefilter('ignore')

    parameter_list = get_parameter_values()
    click.echo('Hyperparameter list {}'.format(parameter_list))

    results = train_elastic_net_model(
        ssgsea_scores_path, phenotypes_path, outer_cv, inner_cv, parameter_list, 'elastic_net', max_iter=max_iterations
    )

    click.echo(results)
