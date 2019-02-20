#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtrellis))


hg38_canonical_seqinfo = readRDS('/repo/201901_locate_discovery_data/annotations/seqinfo_GRCh38.d1.vd1.rds') %>%
    keepStandardChromosomes()

hg38_xlims = tibble(
    chrom = seqnames(hg38_canonical_seqinfo),
    start = 0,  # Same as gtrellis default
    end = seqlengths(hg38_canonical_seqinfo)
) %>%
    filter(chrom != 'chrM')

hg38_cytoband = readRDS('/repo/201901_locate_discovery_data/annotations/cytoband_hg38.rds')


# Cytoband legend
cytoband_legend <- ComplexHeatmap::Legend(
    at = c("p arm", "centromere", "q arm"),
    title = "Cytoband type",
    type = "grid",
    legend_gp = gpar(fille = c("gray50", "red", "gray80")),
    background = 'transparent'
)

# Custom cytoband coloring function to show p/q arms and centromere only
cytoband_coloring <- function(cytoband_chr) {
    cytoband_chr %>%
        transmute(
            col = case_when(
                cytoband_type == 'acen' ~ rgb(217, 47, 39, maxColorValue = 255),
                startsWith(name, 'p') ~ 'gray50',
                startsWith(name, 'q') ~ 'gray80',
                TRUE ~ '#FFFFFF'
            )
        ) %>%
        pull(col)
}

cytoband_panel_drawer <- function(gr) {
    cytoband_chr = gr
    grid.rect(
        cytoband_chr$start, unit(0, "npc"),
        width = cytoband_chr$end - cytoband_chr$start, height = unit(1, "npc"),
        default.units = "native", hjust = 0, vjust = 0,
        gp = gpar(
            fill = cytoband_coloring(cytoband_chr),
            col = 'transparent'
        ))
    # Plot the block
    grid.rect(
        min(cytoband_chr$start), unit(0, "npc"),
        width = max(cytoband_chr$end) - min(cytoband_chr$start), height = unit(1, "npc"),
        default.units = "native", hjust = 0, vjust = 0,
        gp = gpar(fill = "transparent")
    )
}

# Legend of the CNV status
cnv_status_legend <- ComplexHeatmap::Legend(
    at = c("Gain (>0.2)", "Neutral", "Loss (<-0.2)"),
    title = "CN status",
    type = "lines",
    legend_gp = gpar(col = c("red", "gray60", "blue"), lwd = 4),
    background = 'transparent'
)


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
        mutate(
            cnv_status = case_when(
                `log2.copyRatio` >= 0.2 ~ 'Gain',
                `log2.copyRatio` <= -0.2 ~ 'Loss',
                TRUE ~ 'Neutral'
            ),
            # Clamp the min and max value
            `log2.copyRatio` = case_when(
                `log2.copyRatio` > 3 ~ 3,
                `log2.copyRatio` < -3 ~ -3,
                TRUE ~ `log2.copyRatio`
            )
        )
    seg_cnv_gr <- GRanges(
        seqnames = seg_cnv_tbl$chrom,
        ranges = IRanges(start = seg_cnv_tbl$start, end = seg_cnv_tbl$end),
        strand = "*",
        cnv_log2 = seg_cnv_tbl$log2.copyRatio,
        cnv_status = seg_cnv_tbl$cnv_status,
        seqinfo = hg38_canonical_seqinfo
    )
    seg_cnv_gr
}

# Parse command line
args = commandArgs(trailingOnly=TRUE)
sample = args[1]
seg_cnv_pth = args[2]
fig_pth = args[3]


# Read the segment CNV result
seg_cnv_gr = read_seg_cnv(seg_cnv_pth)

pdf(fig_pth, width = 9, height = 6)
# Specify the number of tracks and aesthetics
gtrellis_layout(
    data = hg38_xlims,
    n_track = 3,
    track_axis = c(FALSE, TRUE, FALSE),
    track_height = unit.c(1.5 * grobHeight(textGrob("1XY")), unit(1, "null"), unit(4, 'mm')),
    nrow = 3, compact = TRUE,
    xlab = NULL,
    track_ylab = c('', 'CN log2 ratio', ''),
    track_ylim = c(-3, 3),
    legend = list(cnv_status_legend@grob, cytoband_legend@grob),
    title = str_interp('Somatic copy number ratio of ${sample}')
)
# Track for chromosome names
add_track(panel_fun = function(gr) {
    chr = str_sub(get_cell_meta_data("name"), start = 4)
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr, gp = gpar(fontsize = 10))
})
# Track for CNV ratios
add_segments_track(
    seg_cnv_gr,
    seg_cnv_gr$cnv_log2,
    gp = gpar(
        col = case_when(
            seg_cnv_gr$cnv_status == 'Gain' ~ "red",
            seg_cnv_gr$cnv_status == 'Loss' ~ "blue",
            TRUE ~ "gray60"
        ),
        lwd = 4
    )
)
# Track for cytoband
add_track(
    hg38_cytoband,
    panel_fun = cytoband_panel_drawer
)
dev.off()
