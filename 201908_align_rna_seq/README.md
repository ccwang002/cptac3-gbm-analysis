# 201908_align_rna_seq
Align RNA-seq FASTQs to BAMs following GDC's command


## Installation

    conda create -n gbm_gdc_rna \
        python=3.7 \
        snakemake-minimal=5.5.4 \
        samtools=1.9 htslib=1.9 \
        star=2.6.1a \
        htseq=0.11.2 \
        subread=1.6.4


## Pipeline execution

    # Link all RNA-seq FASTQs
    snakemake link_gdc_rna_fastqs

    # Run STAR alignment
    snakemake star_align_all_samples
    snakemake -j 32 --resources io_heavy=4 -- star_align_all_samples

    # Generate BAM manifests under tracked_results
    snakemake gen_washu_bam_map
