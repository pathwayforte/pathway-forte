# -*- coding: utf-8 -*-

"""This module contain the functional class methods implemented in PathwayForte. For now, GSEA and ssGSEA"""

import itertools as itt

import gseapy
import pandas as pd

from pathway_forte.constants import *
from pathway_forte.mappings import get_mapping_dict, load_compath_mapping_dfs
from pathway_forte.pathway_enrichment.over_representation import log
from pathway_forte.utils import get_num_samples


def create_cls_file(gene_expression_file, normal_sample_file, tumor_sample_file, data):
    """Create categorical (e.g. tumor vs sample) class file format (i.e., .cls) for input into gsea

    :param str gene_expression_file: Text file containing expression values for each gene from each sample.
    :param str normal_sample_file:
    :param str tumor_sample_file:
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

    with open(CLASSES.format(data), 'r') as f:
        data = f.read()
        with open(output, 'w') as w:
            w.write(data[:-1])


def run_gsea(gene_exp, gene_set, phenotype_class, permutations=10, output_dir=GSEA):
    """Run GSEA on a given dataset with a given gene set.

    :param str gene_exp: file with gene expression data
    :param str gene_set: gmt files containing pathway gene sets
    :param str phenotype_class: cls file containing information on class labels
    :param int permutations: number of permutations
    :param str output_dir: output directory
    :return:
    """
    gs_result = gseapy.gsea(
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
        format='png'
    )
    return gs_result


def filter_gsea_results(
        gsea_results_file,
        source,
        kegg_manager=None,
        reactome_manager=None,
        wikipathways_manager=None,
        p_value=None,
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
    if p_value:
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
        wikipathways_manager=wikipathways_manager
    )


def gsea_merge_statistics(merged_pathways_df, dataset):
    """Get statistics for pathways included in the merged gene sets dataFrame. These include the proportion of pathways
    from each of the other databases and the proportion of pathways deriving from 2 or more primary resources

    :param pandas.core.frame.DataFrame merged_pathways_df: dataFrame containing pathways from multiple databases
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

            # add pathway to total number of pathways
            num_of_pathways += 1

            merged_genesets += 1

        if pathway_id.startswith('hsa'):
            num_of_kegg_in_merged += 1
        elif pathway_id.startswith('R-HSA'):
            num_of_reactome_in_merged += 1
        elif pathway_id.startswith('WP'):
            num_of_wikipathways_in_merged += 1
        else:
            raise ValueError('Invalid pathway ID {}.'.format(pathway_id))

    kegg_contributions_to_merged = num_of_kegg_in_merged / num_of_pathways * 100
    reactome_contributions_to_merged = num_of_reactome_in_merged / num_of_pathways * 100
    wikipathways_contributions_to_merged = num_of_wikipathways_in_merged / num_of_pathways * 100
    proportion_of_merged_genesets = merged_genesets / num_of_pathways * 100

    print('For the top {} ranked pathways in the merged {} dataset:'.format(num_of_pathways, dataset))
    print('{0:.2f}% are from KEGG'.format(kegg_contributions_to_merged))
    print('{0:.2f}% are from Reactome'.format(reactome_contributions_to_merged))
    print('{0:.2f}% are from WikiPathways'.format(wikipathways_contributions_to_merged))
    print('{0:.2f}% are a combination of 2 or more databases'.format(proportion_of_merged_genesets))


def rearrange_df_columns(df):
    """Rearrange order of columns"""
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]

    return df


def pathway_names_to_df(filtered_gsea_results_df,
                        all_pathway_ids,
                        source,
                        kegg_manager=None,
                        reactome_manager=None,
                        wikipathways_manager=None
                        ):
    """Get pathway names.

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
        kegg_manager=None,
        reactome_manager=None,
        wikipathways_manager=None,
        p_value=None,
        absolute_nes_filter=None,
        geneset_set_filter_minimum_size=None,
        geneset_set_filter_maximum_size=None
):
    KEGG_GSEA = os.path.join(GSEA, KEGG, 'kegg_{}.tsv').format(dataset)
    REACTOME_GSEA = os.path.join(GSEA, REACTOME, 'reactome_{}.tsv').format(dataset)
    WIKIPATHWAYS_GSEA = os.path.join(GSEA, WIKIPATHWAYS, 'wikipathways_{}.tsv').format(dataset)
    MERGE_GSEA = os.path.join(GSEA, MERGED_GENESET, 'merge_{}.tsv').format(dataset)

    # Load GSEA results and filter dataFrames
    kegg_pathway_df = filter_gsea_results(
        KEGG_GSEA,
        KEGG,
        kegg_manager=kegg_manager,
        p_value=p_value,
        absolute_nes_filter=absolute_nes_filter,
        geneset_set_filter_minimum_size=geneset_set_filter_minimum_size,
        geneset_set_filter_maximum_size=geneset_set_filter_maximum_size
    )
    reactome_pathway_df = filter_gsea_results(
        REACTOME_GSEA,
        REACTOME,
        reactome_manager=reactome_manager,
        p_value=p_value,
        absolute_nes_filter=absolute_nes_filter,
        geneset_set_filter_minimum_size=geneset_set_filter_minimum_size,
        geneset_set_filter_maximum_size=geneset_set_filter_maximum_size

    )
    wikipathways_pathway_df = filter_gsea_results(
        WIKIPATHWAYS_GSEA,
        WIKIPATHWAYS,
        wikipathways_manager=wikipathways_manager,
        p_value=p_value,
        absolute_nes_filter=absolute_nes_filter,
        geneset_set_filter_minimum_size=geneset_set_filter_minimum_size,
        geneset_set_filter_maximum_size=geneset_set_filter_maximum_size
    )
    merged_pathway_df = filter_gsea_results(
        MERGE_GSEA,
        MERGED_GENESET,
        kegg_manager=kegg_manager,
        reactome_manager=reactome_manager,
        wikipathways_manager=wikipathways_manager,
        p_value=p_value,
        absolute_nes_filter=absolute_nes_filter,
        geneset_set_filter_minimum_size=geneset_set_filter_minimum_size,
        geneset_set_filter_maximum_size=geneset_set_filter_maximum_size
    )

    # Merge pathway dataframe without applying filters
    merged_total_df = filter_gsea_results(
        MERGE_GSEA,
        MERGED_GENESET,
        kegg_manager=kegg_manager,
        reactome_manager=reactome_manager,
        wikipathways_manager=wikipathways_manager
    )

    return kegg_pathway_df, reactome_pathway_df, wikipathways_pathway_df, merged_pathway_df, merged_total_df


def get_pairwise_mapping_numbers(
        kegg_pathway_df,
        reactome_pathway_df,
        wikipathways_pathway_df,
        merged_pathway_df,
        merged_total_df
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

    # Get list of pathway IDs from the filtered merged dataFrame
    merged_pathways_actual = merged_pathway_ids_to_list(merged_pathway_df)

    # Get list of pathway IDs from the complete merged dataFrame
    merged_pathways_total = merged_pathway_ids_to_list(merged_total_df)

    # Get pathways in filtered merge dataset with mappings to pathway(s) from other resources
    kegg_pathways_actual = check_pathway_ids(
        kegg_pathway_df, KEGG, merged_pathways_actual, equivalent_mappings_dict
    )
    reactome_pathways_actual = check_pathway_ids(
        reactome_pathway_df, REACTOME, merged_pathways_actual, equivalent_mappings_dict
    )
    wikipathways_pathways_actual = check_pathway_ids(
        wikipathways_pathway_df, WIKIPATHWAYS, merged_pathways_actual, equivalent_mappings_dict
    )

    # Get pathways in total merge dataset with mappings to pathway(s) from other resources
    kegg_pathways_expected = check_pathway_ids(
        kegg_pathway_df, KEGG, merged_pathways_total, equivalent_mappings_dict
    )
    reactome_pathways_expected = check_pathway_ids(
        reactome_pathway_df, REACTOME, merged_pathways_total, equivalent_mappings_dict
    )
    wikipathways_pathways_expected = check_pathway_ids(
        wikipathways_pathway_df, WIKIPATHWAYS, merged_pathways_total, equivalent_mappings_dict
    )

    # Get set of all actual equivalent mappings in filtered merge dataset
    actual_merge = \
        set(tuple(mapping) for mapping in kegg_pathways_actual if mapping) | \
        set(tuple(mapping) for mapping in reactome_pathways_actual if mapping) | \
        set(tuple(mapping) for mapping in wikipathways_pathways_actual if mapping)

    # Get set of all expected equivalent mappings in total merge dataset
    expected_merge = \
        set(tuple(mapping) for mapping in kegg_pathways_expected if mapping) | \
        set(tuple(mapping) for mapping in reactome_pathways_expected if mapping) | \
        set(tuple(mapping) for mapping in wikipathways_pathways_expected if mapping)

    # Get number of actual and expected mappings for merge dataset
    actual_num_dict['Merge dataset'] = len(actual_merge)
    expected_num_dict['Merge dataset'] = len(expected_merge)

    return actual_num_dict, expected_num_dict


def compare_database_results(df_1, resource_1, df_2, resource_2, mapping_dict):
    """Compare pathways in the dataframe from GSEA results to evaluate the concordance in similar pathways."""
    df_1_ids = df_1['pathway_id'].tolist()
    df_2_ids = df_2['pathway_id'].tolist()

    pathways_with_mappings = []
    matching_mappings = []

    for pathway_resource_1 in df_1_ids:

        # Pathway does not have mappings
        if (resource_1, pathway_resource_1) not in mapping_dict:
            continue

        # Get pathway mappings
        mappings_for_pathway = mapping_dict[(resource_1, pathway_resource_1)]

        for resource, mapping_pathway_id in mappings_for_pathway:

            if resource != resource_2:
                continue

            pathways_with_mappings.append(mapping_pathway_id)

            if mapping_pathway_id in df_2_ids:
                matching_mappings.append(mapping_pathway_id)

    return matching_mappings, pathways_with_mappings


def merged_pathway_ids_to_list(df):
    """Split pathway IDs in dataFrame and return list of pathway IDs."""
    return [
        pathway.split('|')
        for pathway in df["pathway_id"].tolist()
    ]


def check_merged_pathways_for_matches(merged_pathways, pathway_id):
    """Check if pathways in the merged data contain a specific pathway ID"""
    for pathway_ids in merged_pathways:
        if pathway_id in pathway_ids:
            if 1 == len(pathway_ids):
                print(pathway_ids)
            return pathway_ids
    return None


def check_pathway_ids(source_df, source, merged_pathways_list, mappings_dict):
    """Check pathway IDs for pathways with mappings"""
    pathways_with_matches = []

    pathway_ids = source_df["pathway_id"].tolist()

    for pathway_id in pathway_ids:
        if (source, pathway_id) not in mappings_dict:
            continue

        pathways_with_matches.append(check_merged_pathways_for_matches(merged_pathways_list, pathway_id))

    return pathways_with_matches


def run_ssgsea(filtered_expression_data, gene_set, output_dir=SSGSEA, processes=1):
    """Run single sample GSEA (ssGSEA) on filtered gene expression data set.

    :param pandas.core.frame.DataFrame expression_data: filtered gene expression values for samples
    :param str gmt_file: .gmt file containing gene sets
    :param output_dir: output directory
    :return: ssGSEA results in respective directory
    """
    ssgsea_result = gseapy.ssgsea(
        data=filtered_expression_data,
        gene_sets=gene_set,
        outdir=output_dir,  # do not write output to disk
        max_size=3000,
        sample_norm_method='rank',  # choose 'custom' for your own rank list
        permutation_num=0,  # skip permutation procedure, because you don't need it
        no_plot=True,  # skip plotting to speed up
        processes=processes, format='png'
    )
    log.info('Done with ssGSEA')
    return ssgsea_result


def filter_gene_exp_data(expression_data, gmt_file):
    """Filter gene expression data file to include only gene names which are found in the gene set files.

    :param pandas.core.frame.DataFrame expression_data: gene expression values for samples
    :param str gmt_file: .gmt file containing gene sets
    :return: Filtered gene expression data with genes with no correspondences in gene sets removed
    :rtype: pandas.core.frame.DataFramekegg_xml_parser.py
    """
    filtered_expression_data = expression_data.copy()

    # Gene universe from gene set
    gene_sets = gseapy.parser.gsea_gmt_parser(gmt_file, max_size=40000)

    gene_universe = set(itt.chain(*gene_sets.values()))

    counter = 0

    genes_to_remove = list()

    for gene in filtered_expression_data.index.values:

        if gene not in gene_universe:
            genes_to_remove.append(gene)
            counter += 1

    log.info(f'Expression data has {len(filtered_expression_data.index.values)}')
    log.info(f'Gene universe has {len(gene_universe)}')
    log.info(f'{counter} were removed in expression data')
    log.info(
        f'{(len(filtered_expression_data.index.values)-counter) /len(gene_universe)*100:.4f}% '
        f'of the gene expression data is mapped to the pathway datasets')

    return filtered_expression_data.drop(genes_to_remove)
