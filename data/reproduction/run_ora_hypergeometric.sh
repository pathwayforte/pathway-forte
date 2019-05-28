#!/usr/bin/env bash

workon pathme

export PF_RESULTS=~/Desktop/pf_results
export PF_REPO=~/dev/pathwayforte/pathway-forte

export INPUT="$PF_RESULTS/input/deseq2"
export GMT_DIR="$PF_REPO/data/gmt_files/analogs"
export OUTPUT="$PF_RESULTS/output/ora"

pathway_forte ora hypergeometric --genesets "$GMT_DIR/kegg_reactome_analogs.gmt" --fold-changes "$INPUT/brca_deseq2.csv" -o "$OUTPUT/ora_geometric_kegg_reactome_brca.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/kegg_reactome_analogs.gmt" --fold-changes "$INPUT/kirc_deseq2.csv" -o "$OUTPUT/ora_geometric_kegg_reactome_kirc.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/kegg_reactome_analogs.gmt" --fold-changes "$INPUT/lihc_deseq2.csv" -o "$OUTPUT/ora_geometric_kegg_reactome_lihc.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/kegg_reactome_analogs.gmt" --fold-changes "$INPUT/prad_deseq2.csv" -o "$OUTPUT/ora_geometric_kegg_reactome_prad.tsv"

pathway_forte ora hypergeometric --genesets "$GMT_DIR/reactome_wikipathways_analogs.gmt" --fold-changes "$INPUT/brca_deseq2.csv" -o "$OUTPUT/ora_geometric_reactome_wikipathways_brca.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/reactome_wikipathways_analogs.gmt" --fold-changes "$INPUT/kirc_deseq2.csv" -o "$OUTPUT/ora_geometric_reactome_wikipathways_kirc.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/reactome_wikipathways_analogs.gmt" --fold-changes "$INPUT/lihc_deseq2.csv" -o "$OUTPUT/ora_geometric_reactome_wikipathways_lihc.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/reactome_wikipathways_analogs.gmt" --fold-changes "$INPUT/prad_deseq2.csv" -o "$OUTPUT/ora_geometric_reactome_wikipathways_prad.tsv"

pathway_forte ora hypergeometric --genesets "$GMT_DIR/wikipathways_kegg_analogs.gmt" --fold-changes "$INPUT/brca_deseq2.csv" -o "$OUTPUT/ora_geometric_wikipathways_kegg_brca.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/wikipathways_kegg_analogs.gmt" --fold-changes "$INPUT/kirc_deseq2.csv" -o "$OUTPUT/ora_geometric_wikipathways_kegg_kirc.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/wikipathways_kegg_analogs.gmt" --fold-changes "$INPUT/lihc_deseq2.csv" -o "$OUTPUT/ora_geometric_wikipathways_kegg_lihc.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/wikipathways_kegg_analogs.gmt" --fold-changes "$INPUT/prad_deseq2.csv" -o "$OUTPUT/ora_geometric_wikipathways_kegg_prad.tsv"

pathway_forte ora hypergeometric --genesets "$GMT_DIR/analogous_mpath.gmt" --fold-changes "$INPUT/brca_deseq2.csv" -o "$OUTPUT/ora_geometric_analogues_mpath_brca.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/analogous_mpath.gmt" --fold-changes "$INPUT/kirc_deseq2.csv" -o "$OUTPUT/ora_geometric_analogues_mpath_kirc.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/analogous_mpath.gmt" --fold-changes "$INPUT/lihc_deseq2.csv" -o "$OUTPUT/ora_geometric_analogues_mpath_lihc.tsv"
pathway_forte ora hypergeometric --genesets "$GMT_DIR/analogous_mpath.gmt" --fold-changes "$INPUT/prad_deseq2.csv" -o "$OUTPUT/ora_geometric_analogues_mpath_prad.tsv"
