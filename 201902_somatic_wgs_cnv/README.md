# 201902_somatic_wgs_cnv
This project plot and analyze the somatic WGS CNVs by the BIC-seq2 output.



## Pipeline execution
Copy the results from the raw pipeline output:

    snakemake copy_bicseq2_cnv_all_samples

Plot the segment CNV under `processed_data/seg_cnv_plot`:

    snakemake plot_seg_cnv_all_samples



## Notebooks
The R Notebooks can be run using the Docker image [lbwang/rocker-genome].

- `plot_segment_cnv.Rmd`: Script to plot the segment CNV of one sample

[lbwang/rocker-transcriptome]: https://hub.docker.com/r/lbwang/rocker-genome
