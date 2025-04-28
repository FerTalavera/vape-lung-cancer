# CNVkit

library(maftools)

gistic_cnvkit = readGistic(gisticDir = "/Users/fernandatalavera/Desktop/BSc_thesis/results/3_copy_number_profile/cnvkit/figures/gistic2/results")

getGeneSummary(gistic_cnvkit)
getSampleSummary(gistic_cnvkit)

CytobandSummary_cnvkit <- getCytobandSummary(gistic_cnvkit)
CytobandSummary_cnvkit <- CytobandSummary_cnvkit[1:10,]
cytoband_ck <- CytobandSummary_cnvkit$Cytoband

amplifications_cnvkit <- read.table("/Users/fernandatalavera/Desktop/BSc_thesis/results/3_copy_number_profile/cnvkit/figures/gistic2/results/amp_genes.conf_90.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
deletions_cnvkit <- read.table("/Users/fernandatalavera/Desktop/BSc_thesis/results/3_copy_number_profile/cnvkit/figures/gistic2/results/del_genes.conf_90.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# CopywriteR

library(maftools)

gistic_copywriter = readGistic(gisticDir = "/Users/fernandatalavera/Desktop/BSc_thesis/results/3_copy_number_profile/copywriter/figures/gistic2/results")

getGeneSummary(gistic_copywriter)
getSampleSummary(gistic_copywriter)

CytobandSummary_copywriter <- getCytobandSummary(gistic_copywriter)
CytobandSummary_copywriter <- CytobandSummary_copywriter[1:10,]
cytobands_cw <- CytobandSummary_copywriter$Cytoband

amplifications_copywriter <- read.table("/Users/fernandatalavera/Desktop/BSc_thesis/results/3_copy_number_profile/copywriter/figures/gistic2/results/amp_genes.conf_90.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
deletions_copywriter <- read.table("/Users/fernandatalavera/Desktop/BSc_thesis/results/3_copy_number_profile/copywriter/figures/gistic2/results/del_genes.conf_90.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Intersections

## cytobands
intersect(cytoband_ck, cytobands_cw)

gisticChromPlot(gistic = gistic_cnvkit, cytobandOffset = 0.06, markBands = c("11qC", "4qD2.2", "4qE2"))
gisticChromPlot(gistic = gistic_copywriter, cytobandOffset = 0.06, markBands = c("11qC", "4qD2.2", "4qE2"))

## genes

### cnvkit
amplifications_selected <- amplifications_cnvkit[, c("X11qC"), drop = FALSE]
deletions_selected <- deletions_cnvkit[, c("X4qD2.2", "X4qE2")]

genes_amp_ck <- amplifications_selected[4:11,1]

genes_del_1 <- deletions_selected[4:6,1]
genes_del_2 <- deletions_selected[4:11,2]
genes_del_ck <- c(genes_del_1, genes_del_2)

### copywriter

amplifications_selected <- amplifications_copywriter[, c("X11qC"), drop = FALSE]
deletions_selected <- deletions_copywriter[, c("X4qD2.2", "X4qE2")]

genes_amp_cw <- amplifications_selected[4:11,1]

genes_del_1 <- deletions_selected[4:6,1]
genes_del_2 <- deletions_selected[4:11,2]
genes_del_cw <- c(genes_del_1, genes_del_2)

intersect(genes_amp_ck,genes_amp_cw)
intersect(genes_del_ck,genes_del_cw)
