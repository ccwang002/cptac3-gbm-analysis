# 201908_gene_quantification
Run gene quantification on GDC and WashU aligned RNA-seq BAMs.


## Installation
Use `gbm_gdc_rna` environment from `201908_align_rna_seq` project.


## Pipeline execution

    # Link all RNA-seq BAMs
    snakemake link_gdc_rna_bams link_washu_rna_bams

    # Run readcounts
    snakemake -j40 --resources io_heavy=10 -- \
        all_featurecounts_stranded_readcount all_htseq_stranded_readcount \
        all_featurecounts_unstranded_readcount all_htseq_unstranded_readcount


## Notebooks (run locally)

