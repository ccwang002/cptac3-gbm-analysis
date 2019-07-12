# 201904_mutation
Collect the somatic mutation  provided by WashU (SomaticWrapper and TinDaisy)
and GDC.


## Installation

    conda create -n gbm_mut \
        'python=3.7' \
        snakemake-minimal sqlalchemy cyvcf2 \
        samtools=1.9 htslib=1.9


## Pipeline execution
Load the somatic mutations by TinDaisy and GDC:

    snakemake link_gdc_raw_vcfs link_gdc_annotated_vcfs
    snakemake all_gdc_mafs all_tindaisy_mafs
    snakemake add_gdc_mutation_to_db add_all_tindaisy_mafs add_all_tindaisy_raw_merged_mafs

