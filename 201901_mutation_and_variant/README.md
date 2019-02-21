# 201901_mutation_and_variant
Collect the somatic mutation and germline variant calls provided by WashU and GDC.


## Installation

    conda create -n gbm_mut \
        'python=3.7' \
        snakemake-minimal sqlalchemy cyvcf2 \
        ipython flake8


## Pipeline execution
Load the somatic mutations by WashU and GDC:

    snakemake add_gdc_mutation_to_db add_washu_mutation_to_db
