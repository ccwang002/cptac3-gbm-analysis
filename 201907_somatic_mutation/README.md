# 201907_somatic_mutation
Collect the somatic mutation generated by WashU (TinDaisy) and GDC.


## Installation

    conda create -n gbm_mut \
        'python=3.7' \
        snakemake-minimal sqlalchemy cyvcf2 \
        samtools=1.9 htslib=1.9

Additionally, download [gdc-client].

[gdc-client]: https://gdc.cancer.gov/access-data/gdc-data-transfer-tool


## Download GDC VCFs

    ./gdc-client download \
        -n 5 --retry-amount 2 \
        -t external_data/gdc_tokens/gdc-user-token.2019-07-15T02_36_07.394Z.txt \
        -m gdc_all_gbm_discovery_vcf_manifest.tsv \
        -d external_data/gdc_client_download


## Pipeline execution
Load the somatic mutations by TinDaisy and GDC:

    # Copy/link raw VCFs
    snakemake --resources ssh_connections=6 -- link_gdc_vcfs link_tindaisy_vcfs

    # Convert and annotation VCFs to be MAFs
    snakemake -j20 all_gdc_mafs all_tindaisy_mafs

    # Combine TinDaisy MAFs into one MAF
    snakemake all_combined_tindaisy_mafs all_combined_gdc_mafs

Once the combined MAFs are generated, remove the individual MAFs by

    rm -rf processed_data/gdc_mafs processed_data/tindaisy_mafs processed_data/tindaisy_raw_annotated_mafs

Currently the database is not used:

    # Load mutations into SQLite database
    snakemake add_gdc_mutation_to_db add_all_tindaisy_mafs add_all_tindaisy_raw_merged_mafs
