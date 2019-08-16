# 201908_align_rna_seq
Align RNA-seq FASTQs to BAMs following GDC's command


## Installation

    conda create -n gbm_gdc_rna \
        python=3.7 \
        snakemake-minimal=5.5.4 \
        samtools=1.9 htslib=1.9 \
        star=2.6.0c


## Pipeline execution

    # Link all RNA-seq FASTQs
    snakemake link_gdc_rna_fastqs

    # Run STAR alignment
    snakemake star_align_all_samples

