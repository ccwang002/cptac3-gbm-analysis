#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))

# Read seqinfo
hg38_seqinfo = readRDS('/repo/annotations/seqinfo_GRCh38.d1.vd1.rds')

# A function to convert the segment TSV of one sample to GRanges
read_seg_cnv <- function(pth) {
    seg_cnv_tbl <- read_tsv(
        pth,
        col_types = cols(
            .default = col_double(),
            chrom = col_character(),
            start = col_integer(),
            end = col_integer(),
            binNum = col_integer()
        )
    ) %>%
        dplyr::select(chrom, start, end, binNum, log2_copy_ratio = log2.copyRatio) %>%
        mutate(
            cnv_status = factor(
                case_when(
                    log2_copy_ratio >= 0.2 ~ 'Gain',
                    log2_copy_ratio <= -0.2 ~ 'Loss',
                    TRUE ~ 'Neutral'
                ),
                levels = c('Loss', 'Neutral', 'Gain')
            )
        )

    # Convert the tibble to GRanges object
    seg_cnv_gr <- makeGRangesFromDataFrame(
        seg_cnv_tbl,
        keep.extra.columns = TRUE,
        ignore.strand = TRUE,
        seqinfo = hg38_seqinfo
    )
    seg_cnv_gr
}

# Parse the command line arguments
args = commandArgs(trailingOnly=TRUE)
# Output GRangesList path
all_seg_cnv_grl_pth <- args[1]
# Path to all samples' segment CNV tsv.gz
seg_cnv_pths <- args[-1]

# Extract the sample names
samples <- str_match(seg_cnv_pths, '(C3[^/]+)\\.tsv\\.gz$')[, 2]

# Merge the per-sample GRanges objects into one GRangesList
all_seg_cnvs_grl <- seg_cnv_pths %>%
    purrr:::map(read_seg_cnv) %>%
    purrr::set_names(nm = samples) %>%
    GRangesList()

saveRDS(all_seg_cnvs_grl, all_seg_cnv_grl_pth)
