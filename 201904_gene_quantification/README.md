# 201904_gene_quantification
RNA-seq gene expression quantification using Salmon.

Three transcript annotations were used:

- GENCODE v30 basic transcripts (CHR region; reference chromosomes only) (Ensembl v96)
- GENCODE v30 comprehensive transcripts (CHR region; reference chromosomes only)
- All Ensembl v96 transcripts

For GTEx normal samples, their original RNA-seq BAMs are copied from MGI first
to extract the reads into paired FASTQ files.


## Installation

    conda create -n gbm_rna \
        python=3.7 \
        snakemake-minimal=5.4.5 \
        samtools=1.9 htslib=1.9 \
        salmon=0.13.1


## Pipeline execution
Generate the Salmon indices of all annotation sources by:

    snakemake build_all_salmon_indices

Link the CPTAC RNA-seq FASTQs from GDC by:

    snakemake link_cptac_gdc_rna_fastqs

To run the gene quantification of all CPTAC samples using all annotation
sources, one will need to change the CPTAC case list `CPTAC_CASE_LIST_PTH` and
`BAM_MAP_PTH` to run different cases on different machines. Then run:

    snakemake salmon_quant_all_cptac_samples

To run the gene quantifications of all GTEx samples, set the `GTEX_FASTQ_FOLDER`:

    snakemake salmon_quant_all_gtex_samples

The gene quantifications will be available at
`processed_data/salmon_quants/{annotation_source}/{sample}/quant.sf`, where
`annotation_source` can be `gencode_basic`, `gencode_comp`, or `ensembl`.


## Notebooks
Notebooks can be run using the Docker image [lbwang/rocker-transcriptome].

[lbwang/rocker-transcriptome]: https://hub.docker.com/r/lbwang/rocker-transcriptome
