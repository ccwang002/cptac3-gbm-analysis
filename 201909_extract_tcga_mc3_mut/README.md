# 201909_extract_tcga_mc3_mut
This project extracts GBM cases only from the GRCh38 lifted TCGA MC3 mutations.


## Pipeline execution
Use pypy3 can speed up the processing, though it only takes less than 3mins to run. Otherwise, use Python 3.5+


    # Use pypy3
    set -x PATH /diskmnt/Projects/Users/lwang/tools/pypy3.6-7.1.1-beta-linux_x86_64-portable/bin $PATH

    pypy3 scripts/extract_maf.py \
        /diskmnt/Datasets/TCGA/MC3/GRCh38_liftOver/mc3.v0.2.8.PUBLIC.GRCh38_converted.maf.gz \
        tcga_gbm_clinical.tsv.gz \
        tcga_tissue_source_site_codes.tsv.gz \
        2> logs/extract_tcga_gbm_mut.log \
        | gzip -9 > processed_data/TCGA_MC3_GRCh38_public_v0.2.8_GBM_only.maf.gz

