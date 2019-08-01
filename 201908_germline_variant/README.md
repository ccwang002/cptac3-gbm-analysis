# 201908_germline_variant
Collect the germline variants.

## Installation
Follows 201907_somatic_mutation's setup


## Pipeline execution

        # Download VCFs and MAFs
        snakemake -j8 --resources ssh_connections=8 -- all_mafs all_vep_vcfs
