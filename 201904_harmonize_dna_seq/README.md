# 201904_harmonize_dna_seq
This project harmonizes the raw sequencing data (WGS and WXS) on GDC.


## Installation

    conda create -n gbm_harmonize \
        python=3.7 snakemake-minimal \
        htslib=1.9 samtools=1.9 \
        biobambam=2.0.87 \
        bwa=0.7.17 \
        picard=2.19.0


## Pipeline execution

### Running snakemake on MGI LSF/bsub queueing system
Use the `./snakemake_via_lsf.sh` bash script to start a bsub job that launch
snakemake in the queue. This ensures the job will continue to run on the cluster.

A example command of running everything is:

    bash snakemake_via_lsf.sh \
        -j10 --resources io_heavy=5 -- \
        all_wxs_picard_mark_dup_bams

where `io_heavy` specifies the maximal concurrent jobs running BAM to fastq,
and `-j10` specifies how many jobs will run at the same time.

### Pipeline steps
Link all the unharmonized WXS BAMs to local folder:

    snakemake link_all_unharmonized_bams

Extract all WXS BAMs to FASTQs:

    snakemake all_wxs_fastqs

Each sample will generate multiple FASTQs under `processed_data/wxs_fqs/` due
to multiple read groups. Note that one should always remove the FASTQs per the sample folder, and the workflow currently doesn't automaticaly clean up this folder.

Align the FASTQs and mark the PCR duplication of all samples by:

    snakemake all_wxs_picard_mark_dup_bams

This step will align all the FASTQs by bwa-mem under
`processed_data/wxs_readgroup_bam`, sort and merge all the readgroup BAMs under
`processed_data/wxs_merged_bam/{sample}.bam`, and finally, mark PCR duplication
by Picard under `processed_data/wxs_mark_dup_bam/{sample}.bam`.

The final output will be at `processed_data/wxs_mark_dup_bam/{sample}.bam`.
