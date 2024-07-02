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

# remove the chr9 centromeres
print('number of rows before removing chr9 centromeres')
print(nrow(ascn))
start_centromere <- 38500001
end_centromere <- 70000000
ascn <- ascn[!(ascn$chr == '9' & ascn$start >= start_centromere & ascn$end <= end_centromere), ]
print('number of rows after removing chr9 centromeres')
print(nrow(ascn))

# subset the clone_ids to only be the intersection of cells that appear in ascn$cell_id
cl <- cl[cl$cell_id %in% ascn$cell_id, ]
print('number of cells in cl after subsetting')
print(length(cl$cell_id))

h1 <- plotHeatmap(ascn,
                  reorderclusters = TRUE,
                  show_clone_label = TRUE,
                  clusters = cl,
                  plottree = FALSE,
                  show_library_label = FALSE,
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


pdf(args$heatmap, width = 11, height = 6)
print(ComplexHeatmap::draw(h1 + h2,
                            ht_gap = unit(0.6, "cm"),
                            column_title = paste0(args$dataset, "\n# of SIGNALS cells: ", length(unique(ascn$cell_id))),
                            column_title_gp = grid::gpar(fontsize = 10),
                            heatmap_legend_side = "bottom",
                            annotation_legend_side = "bottom",
                            show_heatmap_legend = TRUE))
dev.off()
