# 201908_germline_variant
Collect the germline variants.


## Installation
Follows 201907_somatic_mutation's setup


## Pipeline execution

    # Download VCFs and MAFs
    snakemake -j8 --resources ssh_connections=8 -- all_mafs all_vep_vcfs

    # Combine all MAFs
    # Use pypy to speed up the process
    set -x PATH /diskmnt/Projects/Users/lwang/tools/pypy3.6-7.1.1-beta-linux_x86_64-portable/bin $PATH
    pypy3 scripts/combine_mafs.py external_data/germlinewrapper/maf processed_data/combined_mafs/germlinewrapper_all_cases.maf.gz
