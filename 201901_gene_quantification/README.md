# 201901_gene_quantification
RNA-seq gene expression quantification using Salmon.

Two transcript annotations were used:
- GENCODE v29 basic transcripts (Ensembl v94)
- All Ensembl v94 transcripts


## Installation

    conda create -n gbm_rna \
        "python>=3.7" \
        snakemake-minimal \
        salmon ipython

Version of important tools:
- Salmon: 0.12.0


## Pipeline execution
Generate the Salmon indices of all annotation sources by:

    snakemake build_all_salmon_indices

Run the gene quantification of all samples using all annotation sources by:

    snakemake salmon_quant_all_samples
