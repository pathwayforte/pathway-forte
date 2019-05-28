# -*- coding: utf-8 -*-

"""CLI wrapper to import GMT files."""

import logging
import sys
from typing import Iterable, Tuple

from numpy import in1d

__all__ = [
    'gmt_parser',
]

logger = logging.getLogger(__name__)


def spliterate(lines: Iterable[str], sep='\t') -> Iterable[Tuple[str, ...]]:
    """Split each line in the iterable by the given separator."""
    for line in lines:
        yield line.strip().split(sep)


def gmt_parser(
        path: str,
        min_size: int = 3,
        max_size: int = 1000,
        gene_list=None,
        pathway_resource: bool = True,
) -> dict:
    """Parse gene_sets.gmt(gene set database) file or download from enrichr server.

    :param path: the gene_sets.gmt file of GSEA input or an enrichr library name.
     checkout full enrichr library name here: http://amp.pharm.mssm.edu/Enrichr/#stats
    :param min_size: Minimum allowed number of genes from gene set also the data set. Default: 3.
    :param max_size: Maximum allowed number of genes from gene set also the data set. Default: 5000.
    :param gene_list: Used for filtering gene set. Only used this argument for :func:`call` method.
    :return: Return a new filtered gene set database dictionary.
    """
    if not path.lower().endswith(".gmt"):
        raise ValueError("The given file is not a *.gmt file. Please ensure its extension.")

    with open(path) as file:
        if pathway_resource:
            # Special case for PathwayForte since our gene sets include information about the pathway resource too
            genesets_dict = {
                (name, database): genes
                for name, database, *genes in spliterate(file)
            }
            # Normal GMT file: pathway name, description, gene sets...
        else:
            genesets_dict = {
                name: genes
                for name, _, *genes in spliterate(file)
            }

    # Apply gene set filter
    genesets_filter = {
        key: genes
        for key, genes in genesets_dict.items()
        if min_size < len(genes) <= max_size
    }

    if gene_list is not None:
        subsets = sorted(genesets_filter.keys())
        for subset in subsets:
            tag_indicator = in1d(gene_list, genesets_filter.get(subset), assume_unique=True)
            tag_len = sum(tag_indicator)
            if tag_len <= min_size or tag_len >= max_size:
                del genesets_filter[subset]
            else:
                continue
    # some_dict = {key: value for key, value in some_dict.items() if value != value_to_remove}
    # use np.intersect1d() may be faster???
    filsets_num = len(genesets_dict) - len(genesets_filter)
    logging.info(f"{filsets_num} gene sets were removed with filters: max_size={max_size} and min_size={min_size}")

    if filsets_num == len(genesets_dict):
        logging.error("No gene sets passed throught filtering condition!!!, try new paramters again!\n" + \
                      "Note: Gene names for gseapy is case sensitive.")
        sys.exit(1)
    else:
        return genesets_filter
