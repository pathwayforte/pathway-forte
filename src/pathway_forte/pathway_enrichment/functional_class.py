# -*- coding: utf-8 -*-

"""This module contain the functional class methods implemented in PathwayForte. For now, GSEA and ssGSEA"""

import itertools as itt
import logging
import os
from collections import defaultdict
from typing import Optional

import bio2bel_kegg
import bio2bel_reactome
import bio2bel_wikipathways
import gseapy
import pandas as pd
from gseapy.gsea import SingleSampleGSEA
from numpy import sign

from pathway_forte.constants import (
    CLASSES, GSEA, KEGG, MPATH, PHENOTYPE_CLASSES, REACTOME, SSGSEA, WIKIPATHWAYS,
)
from pathway_forte.mappings import get_equivalent_mappings_dict, get_mapping_dict, load_compath_mapping_dfs
from pathway_forte.utils import get_num_samples

logger = logging.getLogger(__name__)


def create_cls_file(gene_expression_file, normal_sample_file, tumor_sample_file, data):
    """Create categorical (e.g. tumor vs sample) class file format (i.e., .cls) for input into GSEA

    :param str gene_expression_file: Text file containing expression values for each gene from each sample.
    :param str normal_sample_file:
    :param str tumor_sample_file:
    :param data:
    """
    num_normal_samples = get_num_samples(normal_sample_file)
    num_tumor_samples = get_num_samples(tumor_sample_file)

    total_num_samples = len(gene_expression_file.columns)
    class_number = 2
    normal_class = 'Normal'
    tumor_class = 'Tumor'

    cls_indices = ['Normal'] * num_normal_samples + ['Tumor'] * num_tumor_samples

    with open(CLASSES.format(data), "w+") as classes_file:
        classes_file.write(
            "{} {} {} \n{} {} {}\n".format(total_num_samples, class_number, 1, "#", normal_class, tumor_class))
        for item in cls_indices:
            classes_file.write("{} ".format(item))

    output = PHENOTYPE_CLASSES.format(data)

    with open(CLASSES.format(data), 'r') as f, open(output, 'w') as w:
        data = f.read()
        w.write(data[:-1])


def run_gsea(gene_exp: str, gene_set: str, phenotype_class: str, permutations: int = 500, output_dir: str = GSEA):
    """Run GSEA on a given dataset with a given gene set.

    :param gene_exp: file with gene expression data
    :param gene_set: gmt files containing pathway gene sets
    :param phenotype_class: cls file containing information on class labels
    :param permutations: number of permutations
    :param output_dir: output directory
    :return:
    """
    return gseapy.gsea(
        data=gene_exp,
        gene_sets=gene_set,
        cls=phenotype_class,  # cls=class_vector
        max_size=3000,
        # set permutation_type to phenotype if samples >=15
        permutation_type='phenotype',
        permutation_num=permutations,  # reduce number to speed up test
        outdir=output_dir,  # do not write output to disk
        no_plot=True,  # Skip plotting
        processes=4,
        format='png',
    )


def filter_gsea_results(
        gsea_results_file,
        source,
        kegg_manager: Optional[bio2bel_kegg.Manager] = None,
        reactome_manager: Optional[bio2bel_reactome.Manager] = None,
        wikipathways_manager: Optional[bio2bel_wikipathways.Manager] = None,
        p_value: Optional[float] = None,
        absolute_nes_filter=None,
        geneset_set_filter_minimum_size=None,
        geneset_set_filter_maximum_size=None,
):
    """Get top and bottom rankings from GSEA results.

    :param gsea_results_file: path to GSEA results in .tsv file format
    :param p_value: maximum p value allowed
    :param absolute_nes_filter: filter by magnitude of normalized enrichment scores
    :param geneset_set_filter_minimum_size: filter to include a minimum number of genes in a gene set
    :param geneset_set_filter_maximum_size: filter to include a maximum number of genes in a gene set
    :param kegg_manager: KEGG manager
    :param reactome_manager: Reactome manager
    :param wikipathways_manager: WikiPathways manager
    :return: list of pathways ranked as having the highest and lowest significant enrichment scores
    :rtype list
    """
    # Read GSEA results to pandas dataFrame
    gsea_results_df = pd.read_csv(gsea_results_file, sep='\t')

    # Filter dataFrame to include only those pathways with a p-value less than X
    if p_value is not None:
        gsea_results_df = gsea_results_df.loc[gsea_results_df['pval'] < p_value]

    #  Filter dataFrame by magnitude of normalized enrichment scores
    if absolute_nes_filter:
        gsea_results_df = gsea_results_df.loc[
            abs(gsea_results_df['nes']) > absolute_nes_filter]

    #  Filter dataFrame by minimun pathway geneset sizes
    if geneset_set_filter_minimum_size:
        gsea_results_df = gsea_results_df.loc[
            gsea_results_df['geneset_size'] > geneset_set_filter_minimum_size]

    #  Filter dataFrame by maximum pathway geneset sizes
    if geneset_set_filter_maximum_size:
        gsea_results_df = gsea_results_df.loc[
            gsea_results_df['geneset_size'] < geneset_set_filter_maximum_size]

    # Sort normalized enrichment scores in descending order
    gsea_results_df = gsea_results_df.sort_values('nes', ascending=False)

    # Rename column
    gsea_results_df = gsea_results_df.rename(columns={'Term': 'pathway_id'})

    # Reset index
    gsea_results_df.reset_index(inplace=True)

    # Remove redundant index column
    del gsea_results_df['index']

    # Get list of pathways after filtering dataFrame
    all_pathway_ids = gsea_results_df['pathway_id'].tolist()

    return pathway_names_to_df(
        gsea_results_df,
        all_pathway_ids,
        source,
        kegg_manager=kegg_manager,
        reactome_manager=reactome_manager,
        wikipathways_manager=wikipathways_manager,
    )


def merge_statistics(merged_pathways_df: pd.DataFrame, dataset: str):
    """Get statistics for pathways included in the merged gene sets dataFrame.

    These include the proportion of pathways from each of the other databases and the proportion of pathways
    deriving from 2 or more primary resources

    :param merged_pathways_df: dataFrame containing pathways from multiple databases
    :return: statistics of contents in merged dataset
    """
    num_of_pathways = len(merged_pathways_df.index)

    # Get the proportion that each database contributes to the merged pathways
    num_of_kegg_in_merged = 0
    num_of_reactome_in_merged = 0
    num_of_wikipathways_in_merged = 0
    merged_genesets = 0

    for pathway_id in merged_pathways_df['pathway_id']:

        if '|' in pathway_id:
            pathway_id.split('|')
            merged_genesets += 1

        if pathway_id.startswith('hsa'):
            num_of_kegg_in_merged += 1
        elif pathway_id.startswith('R-HSA'):
            num_of_reactome_in_merged += 1
        elif pathway_id.startswith('WP'):
            num_of_wikipathways_in_merged += 1
        else:
            raise ValueError(f'Invalid pathway ID {pathway_id}.')

    kegg_contributions_to_merged = num_of_kegg_in_merged / num_of_pathways * 100
    reactome_contributions_to_merged = num_of_reactome_in_merged / num_of_pathways * 100
    wikipathways_contributions_to_merged = num_of_wikipathways_in_merged / num_of_pathways * 100
    proportion_of_merged_genesets = merged_genesets / num_of_pathways * 100

    print('For the {} pathways in the merged dataset results for {}:'.format(num_of_pathways, dataset))
    print('{0:.2f}% are from KEGG'.format(kegg_contributions_to_merged))
    print('{0:.2f}% are from Reactome'.format(reactome_contributions_to_merged))
    print('{0:.2f}% are from WikiPathways'.format(wikipathways_contributions_to_merged))
    print('{0:.2f}% are a combination of 2 or more databases'.format(proportion_of_merged_genesets))


def rearrange_df_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rearrange order of columns"""
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    return df


def get_pathway_names(
        database,
        pathway_df,
        kegg_manager: Optional[bio2bel_kegg.Manager] = None,
        reactome_manager: Optional[bio2bel_reactome.Manager] = None,
        wikipathways_manager: Optional[bio2bel_wikipathways.Manager] = None
):
    if database == KEGG:
        pathway_df['pathway_name'] = [
            kegg_manager.get_pathway_by_id('path:' + pathway_id)
            for pathway_id in list(pathway_df['pathway_id'])
        ]
        return pathway_df

    elif database == REACTOME:
        pathway_df['pathway_name'] = [
            reactome_manager.get_pathway_by_id(pathway_id)
            for pathway_id in list(pathway_df['pathway_id'])
        ]
        return pathway_df

    elif database == WIKIPATHWAYS:
        pathway_df['pathway_name'] = [
            wikipathways_manager.get_pathway_by_id(pathway_id)
            for pathway_id in list(pathway_df['pathway_id'])
        ]
        return pathway_df


def pathway_names_to_df(
        filtered_gsea_results_df,
        all_pathway_ids,
        source,
        kegg_manager: Optional[bio2bel_kegg.Manager] = None,
        reactome_manager: Optional[bio2bel_reactome.Manager] = None,
        wikipathways_manager: Optional[bio2bel_wikipathways.Manager] = None,
) -> pd.DataFrame:
    """Get pathway names.

    :param filtered_gsea_results_df:
    :param all_pathway_ids: list of pathway IDs
    :param source: pathway source (i.e., database name or 'Merged')
    :param kegg_manager: KEGG manager
    :param reactome_manager: Reactome manager
    :param wikipathways_manager: WikiPathways manager
    :return: list of pathway names in rankings
    """
    if source == KEGG:
        filtered_gsea_results_df['pathway_name'] = [
            kegg_manager.get_pathway_by_id('path:' + pathway_id)
            for pathway_id in all_pathway_ids
        ]
        return rearrange_df_columns(filtered_gsea_results_df)

    elif source == REACTOME:
        filtered_gsea_results_df['pathway_name'] = [
            reactome_manager.get_pathway_by_id(pathway_id)
            for pathway_id in all_pathway_ids
        ]
        return rearrange_df_columns(filtered_gsea_results_df)

    elif source == WIKIPATHWAYS:
        filtered_gsea_results_df['pathway_name'] = [
            wikipathways_manager.get_pathway_by_id(pathway_id)
            for pathway_id in all_pathway_ids
        ]
        return rearrange_df_columns(filtered_gsea_results_df)

    merged_rankings = []

    for pathway_ids in all_pathway_ids:
        # A pathway might come from the union of multiple pathways via an equivalent mapping
        equivalent_pathways = []

        for pathway_id in pathway_ids.split('|'):
            if pathway_id.startswith('hsa'):
                kegg_pathway_name = str(kegg_manager.get_pathway_by_id('path:' + pathway_id))
                equivalent_pathways.append(kegg_pathway_name)

            elif pathway_id.startswith('R-HSA'):
                reactome_pathway_name = str(reactome_manager.get_pathway_by_id(pathway_id))
                equivalent_pathways.append(reactome_pathway_name)

            elif pathway_id.startswith('WP'):
                wikipathways_pathway_name = str(wikipathways_manager.get_pathway_by_id(pathway_id))
                equivalent_pathways.append(wikipathways_pathway_name)

        merged_pathways = '|'.join(equivalent_pathways)
        merged_rankings.append(merged_pathways)

    filtered_gsea_results_df['pathway_name'] = merged_rankings

    return rearrange_df_columns(filtered_gsea_results_df)


def gsea_results_to_filtered_df(
        dataset,
        kegg_manager: Optional[bio2bel_kegg.Manager] = None,
        reactome_manager: Optional[bio2bel_reactome.Manager] = None,
        wikipathways_manager: Optional[bio2bel_wikipathways.Manager] = None,
        p_value: Optional[float] = None,
        absolute_nes_filter=None,
        geneset_set_filter_minimum_size=None,
        geneset_set_filter_maximum_size=None
):
    """Get filtered GSEA results dataFrames"""
    kegg_gsea_path = os.path.join(GSEA, KEGG, f'kegg_{dataset}.tsv')
    reactome_gsea_path = os.path.join(GSEA, REACTOME, f'reactome_{dataset}.tsv')
    wikipathways_gsea_path = os.path.join(GSEA, WIKIPATHWAYS, f'wikipathways_{dataset}.tsv')
    merge_gsea_path = os.path.join(GSEA, MPATH, f'merge_{dataset}.tsv')

    # Load GSEA results and filter dataFrames
    kegg_pathway_df = filter_gsea_results(
        kegg_gsea_path,
        KEGG,
        kegg_manager=kegg_manager,
        reactome_manager=reactome_manager,
        wikipathways_manager=wikipathways_manager,
        p_value=p_value,
        absolute_nes_filter=absolute_nes_filter,
        geneset_set_filter_minimum_size=geneset_set_filter_minimum_size,
        geneset_set_filter_maximum_size=geneset_set_filter_maximum_size,
    )
    reactome_pathway_df = filter_gsea_results(
        reactome_gsea_path,
        REACTOME,
        kegg_manager=kegg_manager,
        reactome_manager=reactome_manager,
        wikipathways_manager=wikipathways_manager,
        p_value=p_value,
        absolute_nes_filter=absolute_nes_filter,
        geneset_set_filter_minimum_size=geneset_set_filter_minimum_size,
        geneset_set_filter_maximum_size=geneset_set_filter_maximum_size,

    )
    wikipathways_pathway_df = filter_gsea_results(
        wikipathways_gsea_path,
        WIKIPATHWAYS,
        kegg_manager=kegg_manager,
        reactome_manager=reactome_manager,
        wikipathways_manager=wikipathways_manager,
        p_value=p_value,
        absolute_nes_filter=absolute_nes_filter,
        geneset_set_filter_minimum_size=geneset_set_filter_minimum_size,
        geneset_set_filter_maximum_size=geneset_set_filter_maximum_size,
    )
    merged_pathway_df = filter_gsea_results(
        merge_gsea_path,
        MPATH,
        kegg_manager=kegg_manager,
        reactome_manager=reactome_manager,
        wikipathways_manager=wikipathways_manager,
        p_value=p_value,
        absolute_nes_filter=absolute_nes_filter,
        geneset_set_filter_minimum_size=geneset_set_filter_minimum_size,
        geneset_set_filter_maximum_size=geneset_set_filter_maximum_size,
    )

    # Merge pathway dataframe without applying filters
    merged_total_df = filter_gsea_results(
        merge_gsea_path,
        MPATH,
        kegg_manager=kegg_manager,
        reactome_manager=reactome_manager,
        wikipathways_manager=wikipathways_manager,
        # TODO why not give arguments for other parts?
    )

    return (
        kegg_pathway_df,
        reactome_pathway_df,
        wikipathways_pathway_df,
        merged_pathway_df,
        merged_total_df,
    )


def get_pathways_by_resource(pathways: iter, resource: str) -> list:
    """Return pathways by resource."""
    if resource == KEGG:
        return [
            pathway
            for pathway in pathways
            if pathway.startswith('hsa')
        ]

    elif resource == REACTOME:
        return [
            pathway
            for pathway in pathways
            if pathway.startswith('R-HSA-')
        ]

    elif resource == WIKIPATHWAYS:
        return [
            pathway
            for pathway in pathways
            if pathway.startswith('WP')
        ]

    raise ValueError(f'Pathway database {resource} not found')


def _pairwise_helper(pathways, mappings, source_resource, target_resource):
    """Helper for pairwise comparing pathways."""
    counter = 0

    source_pathways = get_pathways_by_resource(pathways, source_resource)
    target_pathways = get_pathways_by_resource(pathways, target_resource)

    for source_pathway in source_pathways:

        if (source_resource, source_pathway) not in mappings:
            continue

        if any(mapping in target_pathways
               for _, mapping in mappings[(source_resource, source_pathway)]):
            counter += 1

    return counter


def get_analogs_comparison_numbers(
        kegg_reactome_pathway_df,
        reactome_wikipathways_pathway_df,
        wikipathways_kegg_pathway_df,
        *,
        pathway_column="pathway_id"
):
    """Get number of existing versus expected pairwise mappings"""
    # Load mappings
    equivalent_mappings_dict = get_equivalent_mappings_dict()

    actual_num_dict = {}

    # IMPORTANT! These two dictionaries should be keep in the same order as the next one (TO PLOT IN THE RIGHT ORDER)

    actual_num_dict[(KEGG, REACTOME)] = _pairwise_helper(
        kegg_reactome_pathway_df[pathway_column], equivalent_mappings_dict, KEGG, REACTOME
    )

    actual_num_dict[(KEGG, WIKIPATHWAYS)] = _pairwise_helper(
        wikipathways_kegg_pathway_df[pathway_column], equivalent_mappings_dict, KEGG, WIKIPATHWAYS
    )

    actual_num_dict[(REACTOME, KEGG)] = _pairwise_helper(
        kegg_reactome_pathway_df[pathway_column], equivalent_mappings_dict, REACTOME, KEGG
    )

    actual_num_dict[(REACTOME, WIKIPATHWAYS)] = _pairwise_helper(
        reactome_wikipathways_pathway_df[pathway_column], equivalent_mappings_dict, REACTOME, WIKIPATHWAYS
    )

    actual_num_dict[(WIKIPATHWAYS, KEGG)] = _pairwise_helper(
        wikipathways_kegg_pathway_df[pathway_column], equivalent_mappings_dict, WIKIPATHWAYS, KEGG
    )

    actual_num_dict[(WIKIPATHWAYS, REACTOME)] = _pairwise_helper(
        reactome_wikipathways_pathway_df[pathway_column], equivalent_mappings_dict, WIKIPATHWAYS, REACTOME
    )

    # Ensure the number of equivalent pathways are the same
    # assert actual_num_dict[(KEGG, REACTOME)] == actual_num_dict[(REACTOME, KEGG)], \
    #     'Error with KEGG, Reactome'
    #
    # assert actual_num_dict[(KEGG, WIKIPATHWAYS)] == actual_num_dict[(WIKIPATHWAYS, KEGG)], \
    #     'Error with KEGG, WikiPathways'
    #
    # assert actual_num_dict[(WIKIPATHWAYS, REACTOME)] == actual_num_dict[(REACTOME, WIKIPATHWAYS)], \
    #     'Error with Reactome, Wikipathways'

    expected_num_dict = {
        (KEGG, REACTOME): len(
            get_pathways_by_resource(kegg_reactome_pathway_df[pathway_column], KEGG)
        ),
        (KEGG, WIKIPATHWAYS): len(
            get_pathways_by_resource(wikipathways_kegg_pathway_df[pathway_column], KEGG)
        ),
        (REACTOME, KEGG): len(
            get_pathways_by_resource(kegg_reactome_pathway_df[pathway_column], REACTOME)
        ),
        (REACTOME, WIKIPATHWAYS): len(
            get_pathways_by_resource(reactome_wikipathways_pathway_df[pathway_column], REACTOME)
        ),
        (WIKIPATHWAYS, KEGG): len(
            get_pathways_by_resource(wikipathways_kegg_pathway_df[pathway_column], WIKIPATHWAYS)
        ),
        (WIKIPATHWAYS, REACTOME): len(
            get_pathways_by_resource(reactome_wikipathways_pathway_df[pathway_column], WIKIPATHWAYS)
        ),
    }

    return actual_num_dict, expected_num_dict


def get_pairwise_mapping_numbers(
        kegg_pathway_df,
        reactome_pathway_df,
        wikipathways_pathway_df,
):
    """Get number of existing versus expected pairwise mappings"""
    pairwise_comparison = [
        (kegg_pathway_df, KEGG),
        (reactome_pathway_df, REACTOME),
        (wikipathways_pathway_df, WIKIPATHWAYS),
    ]
    # Load mappings
    dfs = load_compath_mapping_dfs()
    final_df = pd.concat([dfs[0], dfs[1], dfs[2]])

    equivalent_mappings_dict = get_mapping_dict(final_df, 'equivalentTo')
    # Dictionary with hierarchical mappings. This dictionary is not used to build the gene set for now, but it could be
    # used in the future for other applications
    part_of_mappings_dict = get_mapping_dict(final_df, 'isPartOf')

    actual_mappings = {}
    expected_mappings = {}

    for (df1, resource1), (df2, resource2) in itt.permutations(pairwise_comparison, 2):
        matching_mappings, pathways_with_mappings = compare_database_results(
            df1, resource1, df2, resource2, equivalent_mappings_dict
        )

        actual_mappings[(resource1, resource2)] = matching_mappings
        expected_mappings[(resource1, resource2)] = pathways_with_mappings

    actual_num_dict = {}
    expected_num_dict = {}

    # Get number of actual mappings between 2 resources
    for resources, mappings in actual_mappings.items():
        actual_num_dict[resources] = len(mappings)

    # Get number of expected mappings between 2 resources
    for resources, mappings in expected_mappings.items():
        expected_num_dict[resources] = len(mappings)

    return actual_num_dict, expected_num_dict


def get_pairwise_mappings(
        kegg_pathway_df,
        reactome_pathway_df,
        wikipathways_pathway_df,
):
    """Get pairwise mappings"""
    pairwise_comparison = [
        (kegg_pathway_df, KEGG),
        (reactome_pathway_df, REACTOME),
        (wikipathways_pathway_df, WIKIPATHWAYS),
    ]
    # Load mappings
    dfs = load_compath_mapping_dfs()
    # Get KEGG-Reactome, KEGG-WikiPathways and WikiPathways-Reactome mappings
    final_df = pd.concat([dfs[0], dfs[1], dfs[2]])

    equivalent_mappings_dict = get_mapping_dict(final_df, 'equivalentTo')

    actual_mappings = {}

    for (df1, resource1), (df2, resource2) in itt.permutations(pairwise_comparison, 2):
        matching_mappings = get_matching_pairs(
            df1, resource1, df2, resource2, equivalent_mappings_dict
        )

        actual_mappings[(resource1, resource2)] = matching_mappings
    return actual_mappings


def compare_database_results(df_1, resource_1, df_2, resource_2, mapping_dict, check_contradiction=False):
    """Compare pathways in the dataframe from enrichment results to evaluate the concordance in similar pathways."""

    # Ensure index is set to pathway id column (not in place)
    df_1 = df_1.set_index('pathway_id')
    df_2 = df_2.set_index('pathway_id')

    # Get pathway ids as lists
    df_1_ids = df_1.index.values
    df_2_ids = df_2.index.values

    pathways_with_mappings = []
    matching_mappings = []

    contradictory_pathways = []

    for pathway_resource_1 in df_1_ids:

        # Pathway does not have mappings
        if (resource_1, pathway_resource_1) not in mapping_dict:
            continue

        # Get all mappings for a pathway
        mappings_for_pathway = mapping_dict[(resource_1, pathway_resource_1)]

        for resource, mapping_pathway_id in mappings_for_pathway:

            if resource != resource_2:
                continue

            # Add all expected pathways to list
            pathways_with_mappings.append(mapping_pathway_id)

            # Add all actual mappings to list
            if mapping_pathway_id in df_2_ids:
                matching_mappings.append(mapping_pathway_id)

                # Check for contradictory results
                if check_contradiction and \
                        sign(df_1.loc[pathway_resource_1]['nes']) != sign(df_2.loc[mapping_pathway_id]['nes']):
                    contradictory_pathways.append((df_1.loc[pathway_resource_1], df_2.loc[mapping_pathway_id]))

    if check_contradiction:
        print(f"Total of #{len(contradictory_pathways)} contradictory pathways: {contradictory_pathways}")

    return matching_mappings, pathways_with_mappings


def get_matching_pairs(df_1, resource_1, df_2, resource_2, equivalent_mappings_dict):
    """Get equivalent pathways and their direction of change."""

    df1_subset = df_1[['pathway_id', 'status']]
    df1_tuples = [tuple(x) for x in df1_subset.values]

    df2_subset = df_2[['pathway_id', 'status']]
    df2_tuples = [tuple(x) for x in df2_subset.values]

    matching_mappings = defaultdict(list)

    for pathway_resource_1, direction_1 in df1_tuples:

        # Pathway does not have mappings
        if (resource_1, pathway_resource_1) not in equivalent_mappings_dict:
            continue

        # Get all mappings for the pathway
        mappings_for_pathway = equivalent_mappings_dict[(resource_1, pathway_resource_1)]

        for resource, mapping_pathway_id in mappings_for_pathway:

            if resource != resource_2:
                continue

            for pathway_resource_2, direction_2 in df2_tuples:

                if pathway_resource_2 == mapping_pathway_id:
                    matching_mappings[resource_1, pathway_resource_1, direction_1].append(
                        (resource_2, pathway_resource_2, direction_2))

    return matching_mappings


def run_ssgsea(
        filtered_expression_data: pd.DataFrame,
        gene_set: str,
        output_dir: str = SSGSEA,
        processes: int = 1,
) -> SingleSampleGSEA:
    """Run single sample GSEA (ssGSEA) on filtered gene expression data set.

    :param filtered_expression_data: filtered gene expression values for samples
    :param gene_set: .gmt file containing gene sets
    :param output_dir: output directory
    :return: ssGSEA results in respective directory
    """
    single_sample_gsea = gseapy.ssgsea(
        data=filtered_expression_data,
        gene_sets=gene_set,
        outdir=output_dir,  # do not write output to disk
        max_size=3000,
        sample_norm_method='rank',  # choose 'custom' for your own rank list
        permutation_num=0,  # skip permutation procedure, because you don't need it
        no_plot=True,  # skip plotting to speed up
        processes=processes,
        format='png',
    )
    logger.info('Done with ssGSEA')
    return single_sample_gsea


def filter_gene_exp_data(expression_data: pd.DataFrame, gmt_file: str):
    """Filter gene expression data file to include only gene names which are found in the gene set files.

    :param expression_data: gene expression values for samples
    :param gmt_file: .gmt file containing gene sets
    :return: Filtered gene expression data with genes with no correspondences in gene sets removed
    :rtype: pandas.core.frame.DataFrame kegg_xml_parser.py
    """
    filtered_expression_data = expression_data.copy()

    # Gene universe from gene set
    gene_sets = gseapy.parser.gsea_gmt_parser(gmt_file, max_size=40000)

    # All the genes in gene set files
    gene_universe = set(itt.chain(*gene_sets.values()))

    genes_to_remove = [
        gene
        for gene in filtered_expression_data.index.values
        if gene not in gene_universe
    ]
    # Genes to be removed because they are not present in the gene sets
    counter = len(genes_to_remove)

    logger.info(f'Expression data has {len(filtered_expression_data.index.values)}')
    logger.info(f'Gene universe has {len(gene_universe)}')
    logger.info(f'{counter} were removed in expression data')
    logger.info(
        f'{(len(filtered_expression_data.index.values) - counter) / len(gene_universe) * 100:.4f}% '
        f'of the gene expression data is mapped to the pathway datasets')

    # Remove non HGNC genes and return dataframe
    return filtered_expression_data.drop(genes_to_remove)
