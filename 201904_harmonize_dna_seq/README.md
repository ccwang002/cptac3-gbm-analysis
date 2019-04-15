# 201904_harmonize_dna_seq
This project harmonizes the raw sequencing data (WGS and WXS) on GDC.


## Installation

    conda create -n gbm_harmonize \
        python=3.7 snakemake-minimal \
        htslib=1.9 samtools=1.9 \
        biobambam=2.0.87 \
        bwa=0.7.17 \
        picard=2.19.0
