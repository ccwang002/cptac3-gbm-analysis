# 201905_somatic_wgs_cnv
This project plot and analyze the somatic WGS CNVs by the BIC-seq2 output.



## Collect the raw output
Copy the raw BIC-seq2 output

    bash scripts/rsync_all_cnv_outputs.sh



## Pipeline execution
Create the R objects of gene and segment CNV of all samples for downstream analysis:

    snakemake merge_all_segment_cnvs merge_all_gene_cnvs

Plot the segment CNV under `processed_data/seg_cnv_plot`:

    snakemake plot_seg_cnv_all_samples



## Notebooks
The R Notebooks can be run using the Docker image [lbwang/rocker-genome].

[lbwang/rocker-genome]: https://hub.docker.com/r/lbwang/rocker-genome
