# -*- coding: utf-8 -*-

"""This module contains all the constants used in the PathwayForte repo."""

import logging
import os
import time

log = logging.getLogger(__name__)

timestr = time.strftime("%d_%m_%Y")

dir_path = os.path.dirname(os.path.realpath(__file__))
SOURCE = os.path.join(os.path.abspath(os.path.join(dir_path, os.pardir)))
DATA = os.path.join(os.path.abspath(os.path.join(SOURCE, os.pardir)), 'data')

"""Cancer Data Sets"""

CANCER_DATA_SETS = {'brca', 'lihc', 'kirc', 'prad', 'ov'}

EXPRESSION_MATRIX = os.path.join(DATA, 'tcga_datasets', '{}', 'expression_matrix_full.txt')
PHENOTYPE_CLASSES = os.path.join(DATA, 'tcga_datasets', '{}', 'phenotype_classes.cls')
CLASSES = os.path.join(DATA, 'tcga_datasets', '{}', 'class.cls')

CLINICAL_DATA = os.path.join(DATA, 'tcga_datasets', '{}', '{}_tcga_clinical_data.tsv')
TUMOR_EXPRESSION_MATRIX = os.path.join(DATA, 'tcga_datasets', '{}', 'tumor_expression_matrix.txt')

NORMAL_EXPRESSION_SAMPLES = os.path.join(DATA, 'tcga_datasets', '{}', 'normal_expression_dimension.txt')
TUMOR_EXPRESSION_SAMPLES = os.path.join(DATA, 'tcga_datasets', '{}', 'tumor_expression_dimension.txt')

RESULTS = os.path.join(DATA, 'results')
CLASSIFIER_RESULTS = os.path.join(RESULTS, 'classifier')


def make_classifier_results_directory():
    """Ensure that the result folder exists."""
    os.makedirs(CLASSIFIER_RESULTS, exist_ok=True)


"""GSEA"""

GSEA = os.path.join(DATA, 'results', 'gsea')

KEGG_GSEA = os.path.join(GSEA, 'kegg')
REACTOME_GSEA = os.path.join(GSEA, 'reactome')
WIKIPATHWAYS_GSEA = os.path.join(GSEA, 'wikipathways')
MERGE_GSEA = os.path.join(GSEA, 'merge')
MSIG_GSEA = os.path.join(GSEA, 'msig')

"""Output files with results for GSEA"""

KEGG_GSEA_TSV = os.path.join(GSEA, 'kegg', 'kegg_{}_{}.tsv')
REACTOME_GSEA_TSV = os.path.join(GSEA, 'reactome', 'reactome_{}_{}.tsv')
WIKIPATHWAYS_GSEA_TSV = os.path.join(GSEA, 'wikipathways', 'wikipathways_{}_{}.tsv')
MERGE_GSEA_TSV = os.path.join(GSEA, 'merge', 'merge_{}_{}.tsv')
CONCATENATED_MERGE_GSEA_TSV = os.path.join(GSEA, 'merge', 'concatenated_merge_{}_{}.tsv')
KEGG_MSIG_GSEA_TSV = os.path.join(GSEA, 'msig', 'msig_kegg_{}_{}.tsv')
REACTOME_MSIG_GSEA_TSV = os.path.join(GSEA, 'msig', 'msig_reactome_{}_{}.tsv')


def make_gsea_export_directories():
    """Ensure that gsea export directories exist."""
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(GSEA, exist_ok=True)
    os.makedirs(KEGG_GSEA, exist_ok=True)
    os.makedirs(REACTOME_GSEA, exist_ok=True)
    os.makedirs(WIKIPATHWAYS_GSEA, exist_ok=True)
    os.makedirs(MERGE_GSEA, exist_ok=True)
    os.makedirs(MSIG_GSEA, exist_ok=True)


"""ssGSEA"""

SSGSEA = os.path.join(DATA, 'results', 'ssgsea')

KEGG_SSGSEA = os.path.join(SSGSEA, 'kegg')
REACTOME_SSGSEA = os.path.join(SSGSEA, 'reactome')
WIKIPATHWAYS_SSGSEA = os.path.join(SSGSEA, 'wikipathways')
MERGE_SSGSEA = os.path.join(SSGSEA, 'merge')
MSIG_SSGSEA = os.path.join(SSGSEA, 'msig')

"""Pickles with results for ssGSEA"""

KEGG_SSGSEA_TSV = os.path.join(SSGSEA, 'kegg', 'kegg_{}_{}.tsv')
REACTOME_SSGSEA_TSV = os.path.join(SSGSEA, 'reactome', 'reactome_{}_{}.tsv')
WIKIPATHWAYS_SSGSEA_TSV = os.path.join(SSGSEA, 'wikipathways', 'wikipathways_{}_{}.tsv')
MERGE_SSGSEA_TSV = os.path.join(SSGSEA, 'merge', 'merge_{}_{}.tsv')
CONCATENATED_MERGE_SSGSEA_TSV = os.path.join(SSGSEA, 'merge', 'concatenated_merge_{}_{}.tsv')
KEGG_MSIG_SSGSEA_TSV = os.path.join(SSGSEA, 'msig', 'kegg_msig_{}_{}.tsv')
REACTOME_MSIG_SSGSEA_TSV = os.path.join(SSGSEA, 'msig', 'reactome_msig_{}_{}.tsv')


def make_ssgsea_export_directories():
    """Ensure that gsea export directories exist."""
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(SSGSEA, exist_ok=True)
    os.makedirs(KEGG_SSGSEA, exist_ok=True)
    os.makedirs(REACTOME_SSGSEA, exist_ok=True)
    os.makedirs(WIKIPATHWAYS_SSGSEA, exist_ok=True)
    os.makedirs(MERGE_SSGSEA, exist_ok=True)
    os.makedirs(MSIG_SSGSEA, exist_ok=True)


"""GMT Files"""

GMT_FOLDER = os.path.join(DATA, 'gmt_files')


def check_gmt_files():
    """Check if GMT files exist and returns GMT files as constant variables."""
    # Get files located in the GMT directory
    gmt_file_names = [
        f for f in os.listdir(GMT_FOLDER)
        if os.path.isfile(os.path.join(GMT_FOLDER, f))
    ]

    # Raise error if files are not found
    if not gmt_file_names:
        log.warning('GMT files missing, please create them by running the "export_gene_sets" command.')
        return None, None, None, None

    kegg_gmt_file, reactome_gmt_file, wikipathways_gmt_file, merge_gmt_file = None, None, None, None

    for file in gmt_file_names:

        if file.startswith('kegg_geneset'):
            kegg_gmt_file = os.path.join(GMT_FOLDER, file)

        elif file.startswith('reactome_geneset'):
            reactome_gmt_file = os.path.join(GMT_FOLDER, file)

        elif file.startswith('wikipathways_geneset'):
            wikipathways_gmt_file = os.path.join(GMT_FOLDER, file)

        elif file.startswith('merged_geneset'):
            merge_gmt_file = os.path.join(GMT_FOLDER, file)

        elif file.startswith('msigdb') or file == 'README.rst' or file == 'concatenated_merge.gmt':
            continue

        else:
            log.warning('Unknown file {} in gmt folder'.format(os.path.join(GMT_FOLDER, file)))

    if not all([kegg_gmt_file, reactome_gmt_file, wikipathways_gmt_file, merge_gmt_file]):
        log.warning('GMT files missing, please create them by running the "export_gene_sets" command.')

    return kegg_gmt_file, reactome_gmt_file, wikipathways_gmt_file, merge_gmt_file


KEGG_GENE_SETS, REACTOME_GENE_SETS, WIKIPATHWAYS_GENE_SETS, MERGED_GENE_SETS = check_gmt_files()

NEW_KEGG_GENE_SETS = os.path.join(GMT_FOLDER, 'kegg_geneset{}.gmt'.format(timestr))
NEW_REACTOME_GENE_SETS = os.path.join(GMT_FOLDER, 'reactome_geneset{}.gmt'.format(timestr))
NEW_WIKIPATHWAYS_GENE_SETS = os.path.join(GMT_FOLDER, 'wikipathways_geneset{}.gmt'.format(timestr))
NEW_MERGED_GENE_SETS = os.path.join(GMT_FOLDER, 'merged_geneset{}.gmt'.format(timestr))

TEMPORAL_KEGG_PATHWAY_GENESET_CSV = os.path.join(GMT_FOLDER, 'kegg_pathway_geneset.csv')
TEMPORAL_REACTOME_PATHWAY_GENESET_CSV = os.path.join(GMT_FOLDER, 'reactome_pathway_geneset.csv')
TEMPORAL_WIKIPATHWAYS_PATHWAY_GENESET_CSV = os.path.join(GMT_FOLDER, 'wikipathways_pathway_geneset.csv')
TEMPORAL_MERGED_PATHWAY_GENESET_CSV = os.path.join(GMT_FOLDER, 'merged_pathway_geneset.csv')

MSIGDB_KEGG_GENE_SETS = os.path.join(GMT_FOLDER, 'msigdb_kegg.gmt')
MSIGDB_REACTOME_GENE_SETS = os.path.join(GMT_FOLDER, 'msigdb_reactome.gmt')
CONCATENATED_MERGE_GENE_SETS = os.path.join(GMT_FOLDER, 'concatenated_merge.gmt')

# Get csv pairwise mapping files
KEGG_REACTOME_URL = "https://raw.githubusercontent.com/ComPath/resources/master/mappings/kegg_reactome.csv"
KEGG_WP_URL = "https://raw.githubusercontent.com/ComPath/resources/master/mappings/kegg_wikipathways.csv"
WP_REACTOME_URL = "https://raw.githubusercontent.com/ComPath/resources/master/mappings/wikipathways_reactome.csv"
SPECIAL_MAPPINGS_URL = "https://raw.githubusercontent.com/ComPath/resources/master/mappings/special_mappings.csv"

RESOURCE = 'Resource'
PATHWAY_ID = 'Pathway ID'
IS_PART_OF = "isPartOf"
MAPPING_TYPE = "Mapping Type"
SOURCE_RESOURCE = "Source Resource"
TARGET_RESOURCE = 'Target Resource'
TARGET_ID = "Target ID"
SOURCE_ID = "Source ID"

KEGG = "kegg"
REACTOME = 'reactome'
WIKIPATHWAYS = 'wikipathways'
MERGED_GENESET = 'merge'
MSIG = 'msig'
CONCATENATED_MERGE = 'concatenated_merge'

PATHWAY_RESOURCES = [KEGG, REACTOME, WIKIPATHWAYS, MERGED_GENESET, MSIG, CONCATENATED_MERGE]

GENESET_COLUMN_NAMES = {
    KEGG: "KEGG Geneset",
    REACTOME: "Reactome Geneset",
    WIKIPATHWAYS: "WikiPathways Geneset",
}

"""Columns to read to perform ORA analysis."""

GENE_SYMBOL = 'gene_symbol'
FOLD_CHANGE = 'log2FoldChange'
P_VALUE = 'p_value'

FC_COLUMNS = {
    GENE_SYMBOL,
    FOLD_CHANGE,
    P_VALUE,
}
