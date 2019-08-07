#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(SummarizedExperiment))


ENSDB_PTH <- '/repo/annotations/EnsDb.Hsapiens.v94.sqlite'
SEQINFO_PTH <- '/repo/annotations/seqinfo_GRCh38.d1.vd1.rds'

# Parse the command line arguments
args = commandArgs(trailingOnly=TRUE)
# Output path
all_gene_cnv_se_pth <- args[1]
# Path to all samples' gene CNV tsv.gz
gene_cnv_pths <- args[-1]


# Read all the samples
samples <- str_match(gene_cnv_pths, '(C3[^/]+)\\.tsv\\.gz$')[, 2]
message(str_interp('Read ${length(samples)} samples.'))

# Read all gene CNV tables and merge into one
gene_cnv_tbl <- samples %>%
    purrr:::map(function(case) {
        pth = str_interp('/repo/201908_somatic_wgs_cnv/processed_data/bicseq2_cnv/gene/${case}.tsv.gz')
        read_tsv(
            pth,
            col_names = c('symbol', 'chrom', 'start', 'end', case),
            col_types = cols(
              symbol = col_character(),
              chrom = col_character(),
              .default = col_double()
            )
        )
    }) %>%
    purrr:::reduce(.f = ~ full_join(.x, .y, by = c('symbol', 'chrom', 'start', 'end')))

# Read seqinfo
hg38_seqinfo <- readRDS(SEQINFO_PTH)

# Read the Ensembl annotation and let it use UCSC chromosome naming style
edb <- EnsDb(ENSDB_PTH)
seqlevelsStyle(edb) <- 'UCSC'

# Select only the protein coding genes on canonical chromosomes.
edb_genes <- genes(
    edb,
    filter = ~ gene_biotype != 'LRG_gene' & seq_name %in% str_c('chr', c(1:22, 'X', 'Y', 'M')),
    columns = c("seq_name", "gene_seq_start", "gene_seq_end", "seq_strand", "gene_id", "gene_id_version", "symbol"),
    return.type = "data.frame"
)

# Try to update the gene symbol by querying the Ensembl annotation
# It will retrieve the Ensembl gene ID which is more stable
new_gene_annotation_tbl <- purrr::pmap_dfr(
    gene_cnv_tbl,
    function(symbol, chrom, start, end, ...) {
        # Try to use the symbol name to find the gene
        g_start = start
        g_end = end
        g_symbol = symbol
        current_gene = edb_genes %>% subset(symbol == g_symbol)
        if (nrow(current_gene) == 1) {
            edb_chrom = current_gene$seq_name
            edb_start = current_gene$gene_seq_start
            edb_end = current_gene$gene_seq_end
            if ((edb_chrom == chrom) & (edb_start == g_start + 1) & (edb_end == g_end)) {
                return(current_gene)
            }
        }
        # Use location to find the gene
        current_gene = edb_genes %>% subset(seq_name == chrom & gene_seq_start == g_start + 1 & gene_seq_end == g_end)
        if (nrow(current_gene) > 1) {
            cat(symbol, 'has more than one entry in EnsDb\n')
            return(NA)
        } else if (nrow(current_gene) == 0) {
            cat(symbol, 'not found\n')
            return(NA)
        }
        return(current_gene)
    }
) %>% as_tibble()

# Add back the original gene symbol per row
row_annotation_tbl <- new_gene_annotation_tbl %>%
    mutate(
        seq_strand = case_when(
            seq_strand == -1 ~ '-',
            seq_strand == 1 ~ '+',
            TRUE ~ '*'
        ),
        original_symbol = gene_cnv_tbl$symbol
    ) %>%
    dplyr::select(-gene_biotype)


# Check if any row has a different gene symbol
rows_with_symbol_change <- row_annotation_tbl %>%
    dplyr::filter(symbol != original_symbol)

message(str_interp('${nrow(rows_with_symbol_change)} rows have different gene symbols after Ensembl query'))

# Create the row gene ranges
row_gr <- makeGRangesFromDataFrame(
    row_annotation_tbl %>% column_to_rownames(var = 'symbol'),
    seqnames.field = 'seq_name',
    start.field = 'gene_seq_start',
    end.field = 'gene_seq_end',
    strand.field = 'seq_strand',
    keep.extra.columns = TRUE,
    seqinfo = hg38_seqinfo
)

# Convert gene CNVs as a matrix
gene_cnv_mat <- as.matrix(gene_cnv_tbl[, -(1:4)])

# Combine all gene CNVs as a SummarizedExperiment object
all_gene_cnv_se <- SummarizedExperiment(
    rowRanges = row_gr,
    assays = list(
        cnv_log2 = gene_cnv_mat
    ),
    metadata = list(
        cohort = 'CPTAC3 GBM adhoc cohort (2019.08)',
        description = 'Gene level somatic CNV (log2) based on whole genome sequencing',
        pipeline = 'BIC-seq2 pipeline at https://github.com/ding-lab/BICSEQ2',
        annotation = 'GENCODE v29 (Ensembl v94) protein coding only',
        contact = 'Liang-Bo Wang <liang-bo.wang@wustl.edu>'
    )
)

saveRDS(all_gene_cnv_se, all_gene_cnv_se_pth)
