#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(gtrellis))


# Read seqinfo of the canonical chromosomes
hg38_canonical_seqinfo = readRDS('/repo/annotations/seqinfo_GRCh38.d1.vd1.rds') %>%
    keepStandardChromosomes()

# Create the genome regions to plot (by default is all the chromosomes without chrM)
hg38_xlims = tibble(
    chrom = seqnames(hg38_canonical_seqinfo),
    start = 0,  # Same as gtrellis default
    end = seqlengths(hg38_canonical_seqinfo)
) %>%
    filter(chrom != 'chrM')

# Read cytoband information
hg38_cytoband = readRDS('/repo/annotations/cytoband_hg38.rds')
# Simplify the cytoband to contain only arms and centromere
hg38_cytoband <- hg38_cytoband %>%
    filter(chrom %in% seqnames(hg38_canonical_seqinfo)) %>%
    mutate(simplified_cytoband_type = case_when(
        cytoband_type == 'acen' ~ 'centromere',
        startsWith(name, 'p') ~ 'p arm',
        startsWith(name, 'q') ~ 'q arm',
        TRUE ~ NA_character_
    )) %>%
    group_by(chrom, simplified_cytoband_type) %>%
    summarize(start = min(start), end = max(end))

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
                simplified_cytoband_type == 'centromere' ~ rgb(217, 47, 39, maxColorValue = 255),
                simplified_cytoband_type == 'p arm' ~ 'gray50',
                simplified_cytoband_type == 'q arm'~ 'gray80',
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

# Parse command line
args = commandArgs(trailingOnly=TRUE)
sample = args[1]
all_seg_cnv_pth = args[2]
fig_pth = args[3]

# Read the segment CNV result
all_seg_cnv_grl = readRDS(all_seg_cnv_pth)

# Select only the sample of interest
seg_cnv_gr = all_seg_cnv_grl[[sample]]

# Clamp the min/max CNV ratio for better visualization
clamp_range = c(min=-2.5, max=2.5)
cn_ratio <- seg_cnv_gr$log2_copy_ratio
cn_ratio[cn_ratio > clamp_range['max']] <- clamp_range['max']
cn_ratio[cn_ratio < clamp_range['min']] <- clamp_range['min']
seg_cnv_gr$log2_copy_ratio <- cn_ratio

# Plot the CNV segments
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
    track_ylim = clamp_range,
    legend = list(cnv_status_legend, cytoband_legend),
    title = str_interp('Somatic copy number ratio of ${sample}')
)
# Track for chromosome names
add_track(panel_fun = function(gr) {
    chr = str_sub(get_cell_meta_data("name"), start = 4)
    grid.rect(gp = gpar(fill = "#EEEEEE"))
    grid.text(chr, gp = gpar(fontsize = 10))
})
# Track for CNV ratios
cnv_color = c("Gain" = 'red', "Loss" = 'blue', "Neutral" = 'gray60')
add_segments_track(
    seg_cnv_gr,
    seg_cnv_gr$log2_copy_ratio,
    gp = gpar(
        col = cnv_color[as.character(seg_cnv_gr$cnv_status)],
        lwd = 4
    )
)
# Track for cytoband
add_track(
    hg38_cytoband,
    panel_fun = cytoband_panel_drawer
)
dev.off()
