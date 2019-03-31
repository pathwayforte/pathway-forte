# -*- coding: utf-8 -*-

"""This module contains the code to export genesets to .gmt files."""

import itertools as itt

import pandas as pd

from pathway_forte.constants import *

__all__ = (
    'get_compath_genesets',
    'create_geneset_df',
    'export_gmt_files',
)


def get_compath_genesets(source, manager):
    """Get ComPath/Bio2BEL gene-sets.

    :param source: either source database or merged
    :param manager: Bio2BEL manager
    :return: dictionary of genesets for pathways from each database: {(source db, pathway id):{pathway_genesets}}
    """
    pathway_genesets = manager.export_gene_sets()

    for key, value in pathway_genesets.items():

        if None in value:
            value.remove(None)

    if source == KEGG:
        return {
            (source, manager.get_pathway_by_name(pathway_name).resource_id.replace('path:', '')): gene_set
            for pathway_name, gene_set in pathway_genesets.items()
            if gene_set
        }

    return {
        (source, manager.get_pathway_by_name(pathway_name).resource_id): gene_set
        for pathway_name, gene_set in pathway_genesets.items()
        if gene_set
    }


def create_geneset_df(all_pathway_genesets, mappings_dict):
    """Create dataframe of genesets for all pathways from all databases. The dataFrame also specifies
     genesets from 2 or more databases for equivalent pathways as merged genesets.

    :param all_pathway_genesets: dictionary of ComPath/Bio2BEL gene sets from all databases
    :param mappings_dict: equivalent pathway pairs dictionary
    :return: dataFrame of all gene sets and their source
    """
    df = pd.DataFrame(columns=GENESET_COLUMN_NAMES.values())

    skip_duplicate_pathways = set()

    # Get genesets for each pathway
    for (database, pathway_id), geneset in all_pathway_genesets.items():

        if (database, pathway_id) in skip_duplicate_pathways:
            continue

        # Geneset for each pathway and its equivalent pathways. If there are no equivalent pathways,
        # the resulting DataFrame will contain only the geneset for the pathway in question.
        row_dict = {PATHWAY_ID: pathway_id, GENESET_COLUMN_NAMES[database]: geneset, RESOURCE: database}

        merged_genesets = [geneset]

        # Check if there are any equivalent mappings for a specific pathway
        if (database, pathway_id) in mappings_dict:

            # Current cell info
            pathway_id_cell = [pathway_id]
            resource_id_cell = [database]

            for (mapping_resource, pathway_mapping_id) in mappings_dict[(database, pathway_id)]:
                equivalent_mapping_genesets = all_pathway_genesets[(mapping_resource, pathway_mapping_id)]

                row_dict[GENESET_COLUMN_NAMES[mapping_resource]] = equivalent_mapping_genesets

                merged_genesets.append(equivalent_mapping_genesets)

                pathway_id_cell.append(pathway_mapping_id)
                resource_id_cell.append(mapping_resource)

                skip_duplicate_pathways.add((mapping_resource, pathway_mapping_id))

            row_dict[PATHWAY_ID] = '|'.join(pathway_id_cell)
            row_dict[RESOURCE] = '|'.join(resource_id_cell)

        merged_genesets_flat = set(itt.chain(*merged_genesets))

        row_dict[MERGED_GENESET] = merged_genesets_flat

        df = df.append(row_dict, ignore_index=True, verify_integrity=True)

    return df


def _remove_equivalence_in_df(df, database):
    """Remove equivalent pathway IDs and resources from dataframe.

    :param df: dataFrame of all gene sets and their source
    :param database: source database or merged if source is combination of pathway resources
    :return: Cleaned dataFrame with duplicate IDs removed
    """
    clean_df = pd.DataFrame(columns=[PATHWAY_ID, RESOURCE, GENESET_COLUMN_NAMES[database]])

    for row_num, row in df.iterrows():

        row_dict = {GENESET_COLUMN_NAMES[database]: row[GENESET_COLUMN_NAMES[database]]}

        resources = row[RESOURCE].split('|')
        pathway_ids = row[PATHWAY_ID].split('|')

        for i in range(0, len(resources)):
            if resources[i] != database:
                continue

            row_dict[RESOURCE] = resources[i]
            row_dict[PATHWAY_ID] = pathway_ids[i]

        if not RESOURCE in row_dict or not PATHWAY_ID in row_dict:
            continue

        clean_df = clean_df.append(row_dict, ignore_index=True, verify_integrity=True)

    return clean_df


def _filter_geneset_file(infile, outfile):
    """Filter gene set file."""
    # Transform .csv to .gmt (Gene Matrix Transposed) file format for gene sets
    special_characters = {
        "{": "",
        "}": "",
        "\'": "",
        ", ": "\t"
    }

    with open(infile, 'r') as f:
        text = f.read()

        for i, j in special_characters.items():
            text = text.replace(i, j)

    with open(outfile, 'w') as w:
        w.write(text)


def export_gmt_files(dataframe):
    """Create dataFrame of genesets for all pathways from all databases. Note that the dataFrame also specifies
    genesets from 2 or more databases for equivalent pathways as merged genesets.
    """
    kegg_df = _remove_equivalence_in_df(dataframe, KEGG)
    reactome_df = _remove_equivalence_in_df(dataframe, REACTOME)
    wikipathways_df = _remove_equivalence_in_df(dataframe, WIKIPATHWAYS)
    merged_geneset_df = dataframe[[PATHWAY_ID, RESOURCE, MERGED_GENESET]]

    kegg_df.to_csv(TEMPORAL_KEGG_PATHWAY_GENESET_CSV, header=False, index=False, sep='\t', encoding='utf-8')
    reactome_df.to_csv(TEMPORAL_REACTOME_PATHWAY_GENESET_CSV, header=False, index=False, sep='\t', encoding='utf-8')
    wikipathways_df.to_csv(
        TEMPORAL_WIKIPATHWAYS_PATHWAY_GENESET_CSV, header=False, index=False, sep='\t', encoding='utf-8'
    )
    merged_geneset_df.to_csv(TEMPORAL_MERGED_PATHWAY_GENESET_CSV, header=False, index=False, sep='\t', encoding='utf-8')

    _filter_geneset_file(TEMPORAL_KEGG_PATHWAY_GENESET_CSV, NEW_KEGG_GENE_SETS)
    _filter_geneset_file(TEMPORAL_REACTOME_PATHWAY_GENESET_CSV, NEW_REACTOME_GENE_SETS)
    _filter_geneset_file(TEMPORAL_WIKIPATHWAYS_PATHWAY_GENESET_CSV, NEW_WIKIPATHWAYS_GENE_SETS)
    _filter_geneset_file(TEMPORAL_MERGED_PATHWAY_GENESET_CSV, NEW_MERGED_GENE_SETS)

    # Remove intermediate csv file as .gmt file contains all genesets and is sole file needed
    os.remove(TEMPORAL_KEGG_PATHWAY_GENESET_CSV)
    os.remove(TEMPORAL_REACTOME_PATHWAY_GENESET_CSV)
    os.remove(TEMPORAL_WIKIPATHWAYS_PATHWAY_GENESET_CSV)
    os.remove(TEMPORAL_MERGED_PATHWAY_GENESET_CSV)
