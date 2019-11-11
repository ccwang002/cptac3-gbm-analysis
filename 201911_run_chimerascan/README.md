# 201911_run_chimerascan
Call RNA gene fusion using chimerascan.

## Installation

    docker pull cgrlab/chimerascan

## Build Index

    docker run -it --rm -u (id -u):(id -g) \
        --group-add (grep cptac /etc/group | cut -d: -f3) \
        -v /diskmnt:/diskmnt \
        -v $PWD:/pwd \
        cgrlab/chimerascan \
        chimerascan_index.py \
            /diskmnt/Datasets/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa \
            /pwd/external_data/gdc_gencode_v22_flatten.tsv \
            /diskmnt/Projects/cptac_scratch_4/CPTAC3_GBM/201911_chimerascan/GDC_GENCODE_v22_index