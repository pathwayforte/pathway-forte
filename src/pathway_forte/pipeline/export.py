# -*- coding: utf-8 -*-

"""CLI wrapper to export updated gene sets using ComPath."""
import logging

from pathway_forte.export_genesets_to_gmt import create_geneset_df, export_gmt_files, get_all_pathway_genesets
from pathway_forte.mappings import get_equivalent_mappings_dict

__all__ = [
    'do_export',
]

logger = logging.getLogger(__name__)


def do_export():
    all_pathway_genesets = get_all_pathway_genesets()
    equivalent_mappings_dict = get_equivalent_mappings_dict()
    df = create_geneset_df(all_pathway_genesets, equivalent_mappings_dict)
    export_gmt_files(df)
    logger.info('Done creating files')
