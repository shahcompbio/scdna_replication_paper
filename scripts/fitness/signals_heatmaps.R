library(argparse)
library(data.table)
library(tidyverse)
library(signals)
library(ggplot2)


parser <- ArgumentParser()
parser$add_argument("--ascn", default=NULL, type="character",
                    help="csv file of signals results")
parser$add_argument("--clones", default=NULL, type="character",
                    help="tsv file that matches cell and clone IDs")
parser$add_argument("--dataset", default=NULL, type="character",
                    help="dataset ID")
parser$add_argument("--heatmap", default=NULL, type="character",
                    help="heatmap plot")
args <- parser$parse_args()
print(args)


# load the ascn data
ascn <- fread(args$ascn)
print(ascn)

# load the clone_ids for each cell
cl <- fread(args$clones, sep='\t')
print(cl)

h1 <- plotHeatmap(ascn,
                  reorderclusters = TRUE,
                  show_clone_label = TRUE,
                  clusters = cl,
                  plottree = FALSE,
                  plotcol = "state")
print('Made h1')
h2 <- plotHeatmap(ascn,
                  reorderclusters = TRUE,
                  show_clone_label = FALSE,
                  clusters = cl,
                  plottree = FALSE,
                  show_library_label = FALSE,
                  plotcol = "state_phase")
print('Made h2')


png(args$heatmap, width = 25, height = 12, units = "in", res = 300)
print(ComplexHeatmap::draw(h1 + h2,
                            ht_gap = unit(0.6, "cm"),
                            column_title = paste0(args$dataset, "\nTotal number of cells: ", length(unique(ascn$cell_id))),
                            column_title_gp = grid::gpar(fontsize = 20),
                            heatmap_legend_side = "bottom",
                            annotation_legend_side = "bottom",
                            show_heatmap_legend = TRUE))
dev.off()
