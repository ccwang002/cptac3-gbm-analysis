# 201903_methylation_array
This project post-processes the Illumina EPIC DNA methylation microarray data from Wen-wei and Sunantha's DNA methylation pipeline.


## Installation

    conda create -n gbm_meth \
        python=3.7 \
        pandas feather-format


## CPTAC Methylation Data (batch 1 to 9)

    rsync -a \
        /diskmnt/Projects/Users/wliang/Methylome_CPTAC/01_Data/GBM/GBM_allPatients.feather \
        processed_data/cptac3_gbm_epic_meth_array.batch1-9.feather


## Annotations
It requires the EPIC probe manifest and annotation in hg38 obtained from <https://zwdzwd.github.io/InfiniumAnnotation> (DOI: 10.1093/nar/gkw967).

    curl --create-dirs -o external_data/EPIC.hg38.manifest.rds -L \
        http://zwdzwd.io/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.rds
    curl --create-dirs -o external_data/EPIC.hg38.manifest.gencode.v22.rds -L \
        http://zwdzwd.io/InfiniumAnnotation/20180909/EPIC/EPIC.hg38.manifest.gencode.v22.rds


## Notebooks
The R Notebooks can be run using the Docker image [lbwang/rocker-genome].

[lbwang/rocker-genome]: https://hub.docker.com/r/lbwang/rocker-genome
