# 201908_germline_variant
Collect the germline variants.


## Installation
Follows 201907_somatic_mutation's setup


## Pipeline execution

    # Download VCFs and MAFs
    mkdir -p external_data/germlinewrapper/vep_vcf
    rsync -a --info=progress2 \
        vw5.gsc.wustl.edu:'/gscmnt/gc2508/dinglab/nvterekhanova/cptac3/germline_calling/prepare_files/Scripts/fernanda/VEP-annotation/out/{57_GBMs,erik_gbm_hg38}/*/*.vep95.vcf' \
        external_data/germlinewrapper/vep_vcf/
    parallel --bar -j8 'bgzip -l9 {}' ::: external_data/germlinewrapper/vep_vcf/*.vep95.vcf

    # Convert VCF to MAF
    snakemake -j 10 all_vcfs all_mafs
