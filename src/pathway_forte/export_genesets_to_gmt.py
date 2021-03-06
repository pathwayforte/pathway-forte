# -*- coding: utf-8 -*-

"""This module contains the code to export genesets to .gmt files."""

import itertools as itt
from _collections import defaultdict

import bio2bel_kegg
import bio2bel_reactome
import bio2bel_wikipathways
import pandas as pd
from compath_utils import CompathManager
from pathway_forte.constants import *
from pathway_forte.utils import handle_file_download, handle_zipfile_download

__all__ = [
    'get_all_pathway_genesets',
    'get_compath_genesets',
    'create_geneset_df',
    'export_gmt_files',
]

logger = logging.getLogger(__name__)


def get_all_pathway_genesets():
    """Get ComPath/Bio2BEL gene sets for each pathway, from each database."""
    logger.info('Getting ComPath/Bio2BEL KEGG gene sets')
    kegg_manager = bio2bel_kegg.Manager()
    kegg_gene_set = get_compath_genesets(KEGG, kegg_manager)

    logger.info('Getting ComPath/Bio2BEL Reactome gene sets')
    reactome_manager = bio2bel_reactome.Manager()
    reactome_gene_set = get_compath_genesets(REACTOME, reactome_manager)

    logger.info('Getting ComPath/Bio2BEL WikiPathways gene sets')
    wikipathways_manager = bio2bel_wikipathways.Manager()
    wikipathways_gene_set = get_compath_genesets(WIKIPATHWAYS, wikipathways_manager)

    # Get ComPath/Bio2BEL genesets for each pathway, from each database
    return {
        **kegg_gene_set,
        **reactome_gene_set,
        **wikipathways_gene_set,
    }


def get_compath_genesets(source, manager: CompathManager):
    """Get ComPath/Bio2BEL gene-sets.

    :param source: either source database or merged
    :param manager: Bio2BEL manager
    :return: dictionary of genesets for pathways from each database: {(source db, pathway id):{pathway_genesets}}
    """
    pathway_genesets = manager.export_gene_sets()

    # Remove empty pathways (pathways without any gene associated to it)
    for key, value in pathway_genesets.items():
        if None in value:
            value.remove(None)

    # Special case for KEGG (prefix needs to be removed)
    if source == KEGG:
        return {
            (source, manager.get_pathway_by_name(pathway_name).resource_id.replace('path:', '')): gene_set
            for pathway_name, gene_set in pathway_genesets.items()
            if gene_set
        }

    # Convert dictionary from Pathway ID to Pathway Name
    return {
        (source, manager.get_pathway_by_name(pathway_name).resource_id): gene_set
        for pathway_name, gene_set in pathway_genesets.items()
        if gene_set
    }


def create_geneset_df(all_pathway_genesets, mappings_dict):
    """Create dataframe of gene sets for all pathways from all databases using the mappings.

    The dataFrame also specifies genesets from 2 or more databases for equivalent pathways as merged genesets.

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
        row_dict = {
            PATHWAY_ID: pathway_id,
            GENESET_COLUMN_NAMES[database]: geneset,
            RESOURCE: database,
        }

        merged_genesets = [geneset]

        # Check if there are any equivalent mappings for a specific pathway
        if (database, pathway_id) in mappings_dict:
            # Current cell info
            pathway_id_cell = [pathway_id]
            resource_id_cell = [database]

            for (mapping_resource, pathway_mapping_id) in mappings_dict[(database, pathway_id)]:
                equivalent_mapping_genesets = all_pathway_genesets[mapping_resource, pathway_mapping_id]

                row_dict[GENESET_COLUMN_NAMES[mapping_resource]] = equivalent_mapping_genesets

                merged_genesets.append(equivalent_mapping_genesets)

                pathway_id_cell.append(pathway_mapping_id)
                resource_id_cell.append(mapping_resource)

                skip_duplicate_pathways.add((mapping_resource, pathway_mapping_id))

            row_dict[PATHWAY_ID] = '|'.join(pathway_id_cell)
            row_dict[RESOURCE] = '|'.join(resource_id_cell)

        merged_genesets_flat = set(itt.chain(*merged_genesets))

        row_dict[MPATH] = merged_genesets_flat

        df = df.append(row_dict, ignore_index=True, verify_integrity=True)

    return df


def _remove_equivalence_in_df(df: pd.DataFrame, database: str) -> pd.DataFrame:
    """Remove equivalent pathway IDs and resources from dataframe.

    :param df: dataFrame of all gene sets and their source
    :param database: source database or merged if source is combination of pathway resources
    :return: Cleaned dataFrame with duplicate IDs removed
    """
    clean_df = pd.DataFrame(columns=[PATHWAY_ID, RESOURCE, GENESET_COLUMN_NAMES[database]])

    for row_num, row in df.iterrows():
        row_dict = {
            GENESET_COLUMN_NAMES[database]: row[GENESET_COLUMN_NAMES[database]]
        }

        resources = row[RESOURCE].split('|')
        pathway_ids = row[PATHWAY_ID].split('|')

        for resource, pathway_id in zip(resources, pathway_ids):
            if resource != database:
                continue

            row_dict[RESOURCE] = resource
            row_dict[PATHWAY_ID] = pathway_id

        if RESOURCE not in row_dict or PATHWAY_ID not in row_dict:
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
        ", ": "\t",
    }

    with open(infile, 'r') as file:
        text = file.read()

    for old, new in special_characters.items():
        text = text.replace(old, new)

    with open(outfile, 'w') as file:
        file.write(text)


def export_gmt_files(df: pd.DataFrame):
    """Create dataFrame of gene sets for all pathways from all databases.

    Note that the dataframe also specifies gene sets from 2 or more databases for equivalent pathways as merged gene
    sets.
    """
    kegg_df = _remove_equivalence_in_df(df, KEGG)
    reactome_df = _remove_equivalence_in_df(df, REACTOME)
    wikipathways_df = _remove_equivalence_in_df(df, WIKIPATHWAYS)
    merged_geneset_df = df[[PATHWAY_ID, RESOURCE, MPATH]]

    # TODO replace these temp paths with proper usage of tempfile
    kegg_df.to_csv(TEMP_KEGG_PATHWAY_GENESET_CSV, header=False, index=False, sep='\t', encoding='utf-8')
    reactome_df.to_csv(TEMP_REACTOME_PATHWAY_GENESET_CSV, header=False, index=False, sep='\t', encoding='utf-8')
    wikipathways_df.to_csv(TEMP_WIKIPATHWAYS_PATHWAY_GENESET_CSV, header=False, index=False, sep='\t', encoding='utf-8')
    merged_geneset_df.to_csv(TEMP_MERGED_PATHWAY_GENESET_CSV, header=False, index=False, sep='\t', encoding='utf-8')

    _filter_geneset_file(TEMP_KEGG_PATHWAY_GENESET_CSV, NEW_KEGG_GENE_SETS)
    _filter_geneset_file(TEMP_REACTOME_PATHWAY_GENESET_CSV, NEW_REACTOME_GENE_SETS)
    _filter_geneset_file(TEMP_WIKIPATHWAYS_PATHWAY_GENESET_CSV, NEW_WIKIPATHWAYS_GENE_SETS)
    _filter_geneset_file(TEMP_MERGED_PATHWAY_GENESET_CSV, NEW_MERGED_GENE_SETS)

    # Remove intermediate csv file as .gmt file contains all genesets and is sole file needed
    os.remove(TEMP_KEGG_PATHWAY_GENESET_CSV)
    os.remove(TEMP_REACTOME_PATHWAY_GENESET_CSV)
    os.remove(TEMP_WIKIPATHWAYS_PATHWAY_GENESET_CSV)
    os.remove(TEMP_MERGED_PATHWAY_GENESET_CSV)


def export_pathbank_geneset(outfile: str):
    """Download PathBank content and export PathBank gene sets for pathways that have mappings to PathMe."""

    logger.info('Downloading PathBank content.')
    # Get PathBank pathways
    handle_zipfile_download(PATHBANK_PATHWAYS_PATH)

    # Get PathBank proteins
    handle_zipfile_download(PATHBANK_PROTEINS_PATH)

    logger.info('Getting ComPath mappings.')
    # Get PathBank-PathMe mappings
    handle_file_download(PATHBANK_KEGG_MAPPINGS, PATHBANK_KEGG_FILE)
    handle_file_download(PATHBANK_REACTOME_MAPPINGS, PATHBANK_REACTOME_FILE)
    handle_file_download(PATHBANK_WIKIPATHWAYS_MAPPINGS, PATHBANK_WIKIPATHWAYS_FILE)

    # Read pathway ID mapping file
    df_pathway_ids = pd.read_csv(PATHBANK_PATHWAYS_FILE, sep=',')
    df_pathway_ids = df_pathway_ids.rename(columns={'SMPDB ID': 'SMPDBID', 'PW ID': 'PWID'})

    # Read pathway-gene membership file
    df_pathway_genes = pd.read_csv(PATHBANK_PROTEINS_FILE, sep=',')

    # Get relevant columns
    df_id_name_gene = df_pathway_genes[['PathBank ID', 'Pathway Name', 'Gene Name']]
    df_id_name_gene = df_id_name_gene.rename(
        columns={
            'PathBank ID': 'PathBankID',
            'Pathway Name': 'PathwayName',
            'Gene Name': 'GeneName'
        }
    )

    # Drop genes if they are NaN
    df_id_name_gene.dropna(subset=['GeneName'], inplace=True)

    df_pathways_merged = pd.merge(
        df_id_name_gene,
        df_pathway_ids,
        left_on='PathBankID',
        right_on='SMPDBID',
    )

    df_pathways_merged = df_pathways_merged[['PWID', 'PathwayName', 'GeneName']]

    # Concatenate PathBank mapping files
    pathbank_mapping_files = [PATHBANK_KEGG_FILE, PATHBANK_REACTOME_FILE, PATHBANK_WIKIPATHWAYS_FILE]
    with open('pathbank_mappings.csv', 'w') as f:
        for fname in pathbank_mapping_files:
            with open(fname) as infile:
                f.write(infile.read())

    pathbank_mappings = {}

    # Get PathBank pathways with mappings to PathMe
    all_pathbank_mappings_df = pd.read_csv('pathbank_mappings.csv', sep=',')
    all_pathbank_mappings_df.rename(
        columns={
            'Source Resource': 'SourceResource',
            'Source ID': 'SourceID',
            'Source Name': 'SourceName',
            'Mapping Type': 'MappingType',
            'Target Resource': 'TargetResource',
            'Target ID': 'TargetID',
            'Target Name': 'TargetName'
        },
        inplace=True,
    )

    pathbank_mappings = {}

    for row in all_pathbank_mappings_df.itertuples(index=False):
        if row.SourceResource == 'pathbank':
            pathbank_mappings[row.SourceName] = row.SourceID
        elif row.TargetResource == 'pathbank':
            pathbank_mappings[row.TargetName] = row.TargetID

    geneset_dict = defaultdict(set)

    # Get pathway-gene set dictionary
    for row in df_pathways_merged.itertuples(index=False):

        for name, identifier in pathbank_mappings.items():

            if row.PWID == identifier:
                geneset_dict[identifier].add(row.GeneName)

    geneset_df = pd.DataFrame.from_dict(data=geneset_dict, orient='index')

    # Add resource column
    geneset_df['Resource'] = pd.Series('pathbank', index=geneset_df.index)
    geneset_df = geneset_df[['Resource'] + [col for col in geneset_df.columns if col != 'Resource']]

    logger.info('Exporting PathBank gene sets.')

    # Write PathBank genesets to gmt file format
    geneset_df.to_csv('pathbank_file.csv', header=False, sep='\t')

    with open('pathbank_file.csv', 'r') as file:
        with open(outfile, 'w') as f:
            for line in file:
                line = line.rstrip()
                f.write(line + '\n')
    f.close()
