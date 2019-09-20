# 201909_call_structural_variant
Call structural variants (SVs) 


## Installation

    conda create -n gbm_sv \
        python=3.7 \
        snakemake-minimal=5.6 \
        samtools=1.9 htslib=1.9 bcftools=1.9

Additional packages of manta are managed by snakemake.

Install the static binary of DELLY v0.8.1 that has OpenMP support:

    curl -Lo delly https://github.com/dellytools/delly/releases/download/v0.8.1/delly_v0.8.1_linux_x86_64bit
    chmod +x delly

And download the regions to exclude SV calling from DELLY code repository.

    curl -Lo external_data/delly/human.hg38.excl.tsv \
        https://raw.githubusercontent.com/dellytools/delly/master/excludeTemplates/human.hg38.excl.tsv
