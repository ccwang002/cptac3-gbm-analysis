# 201905_somatic_wgs_cnv
This project plot and analyze the somatic WGS CNVs by the BIC-seq2 output.



## Collect the raw output
Copy the raw BIC-seq2 output

    rsync -az --info=progress2 \
        --exclude-from=cnv_rsync_exclude_patterns.list \
        vw3.gsc.wustl.edu:/gscmnt/gc2541/cptac3_analysis/gbm_somatic_wgs_cnv/results.GBM-custom-harm.MGI/ \
        external_data/bicseq2_cnv

    rsync -a --info=progress2 \
        --exclude-from=cnv_rsync_exclude_patterns.list \
        /diskmnt/Projects/cptac_downloads_4/CPTAC3_GBM/201905_somatic_wgs_cnv/results.GBM-not-in-discover.katmai/ \
        external_data/bicseq2_cnv

    rsync -a --info=progress2 \
        --exclude-from=cnv_rsync_exclude_patterns.list \
        /diskmnt/Projects/cptac_downloads_4/CPTAC3_GBM/201902_somatic_wgs_cnv/run_cases.GBM-subset/ \
        external_data/bicseq2_cnv

Compress the output

    parallel --bar -j8 'gzip -9 {}' \
        ::: external_data/bicseq2_cnv/*/segmentation/*.cnv \
        external_data/bicseq2_cnv/*/annotation/*.gene_level.log2.seg



## Pipeline execution
Create the R objects of gene and segment CNV of all samples for downstream analysis:

    snakemake merge_all_segment_cnvs merge_all_gene_cnvs

Plot the segment CNV under `processed_data/seg_cnv_plot`:

    snakemake plot_seg_cnv_all_samples



## Notebooks
The R Notebooks can be run using the Docker image [lbwang/rocker-genome].

[lbwang/rocker-genome]: https://hub.docker.com/r/lbwang/rocker-genome
