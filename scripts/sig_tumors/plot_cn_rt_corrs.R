library(argparse)
library(data.table)
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)


parser <- ArgumentParser()
parser$add_argument("--sample_rt_corrs", default=NULL, type="character",
                    help="csv file of sample RT correlations")
parser$add_argument("--sample_cn_dists", default=NULL, type="character",
                    help="csv file of sample CN distances")
parser$add_argument("--clone_rt_corrs", default=NULL, type="character",
                    help="csv file of clone RT correlations")
parser$add_argument("--clone_cn_dists", default=NULL, type="character",
                    help="csv file of clone CN distances")
parser$add_argument("--datasets", default=NULL, type="character", nargs='+',
                    help="list of datasets")
parser$add_argument("--types", default=NULL, type="character", nargs='+',
                    help="list of cancer (tissue) types for each dataset")
parser$add_argument("--signatures", default=NULL, type="character", nargs='+',
                    help="list of signatures for each dataset")
parser$add_argument("--sample_heatmap", default=NULL, type="character",
                    help="sample heatmap plot")
parser$add_argument("--clone_heatmap", default=NULL, type="character",
                    help="clone heatmap plot")
args <- parser$parse_args()
print(args)


# load the correlation and distance matrices
sample_rt_corrs <- fread(args$sample_rt_corrs)
sample_cn_dists <- fread(args$sample_cn_dists)
clone_rt_corrs <- fread(args$clone_rt_corrs)
clone_cn_dists <- fread(args$clone_cn_dists)
print(sample_rt_corrs)

# create a data frame with the datasets, cancer types, and signatures
datasets_df <- data.frame(dataset = args$datasets,
                       type = args$types,
                       signature = args$signatures)

# extract the order that each sample appears in the correlation and distance matrices
# it is the first prefix of the V1 column
sample_rt_corrs$dataset <- gsub("_.*", "", sample_rt_corrs$V1)
sample_cn_dists$dataset <- gsub("_.*", "", sample_cn_dists$V1)
clone_rt_corrs$dataset <- gsub("_.*", "", clone_rt_corrs$V1)
clone_cn_dists$dataset <- gsub("_.*", "", clone_cn_dists$V1)

print(sample_rt_corrs$dataset)
print(sample_cn_dists$dataset)

# check that the order of the samples in the correlation and distance matrices match
stopifnot(all(sample_rt_corrs$dataset == sample_cn_dists$dataset))
stopifnot(all(clone_rt_corrs$dataset == clone_cn_dists$dataset))

# sort the rows of datasets_df to match the order of the rows in the correlation and distance matrices
print(datasets_df)
datasets_df <- datasets_df[match(sample_rt_corrs$dataset, datasets_df$dataset),]
print(datasets_df)

# remove the dataset and V1 columns from the correlation and distance matrices
sample_rt_corrs <- sample_rt_corrs[, -c('dataset', 'V1')]
sample_cn_dists <- sample_cn_dists[, -c('dataset', 'V1')]

# merge datasets_df with the clone correlation and distance matrices
# this will add the cancer type and signature to each clone
clone_rt_corrs <- merge(clone_rt_corrs, datasets_df, by = "dataset")
clone_cn_dists <- merge(clone_cn_dists, datasets_df, by = "dataset")

# extract a clone_df which only has the dataset, cancer type, and signature for each row in the clone correlation and distance matrices
clone_df <- clone_rt_corrs[, c('dataset', 'type', 'signature')]
print(clone_df)

# remove the dataset and V1 columns from the clone correlation and distance matrices
print(colnames(clone_rt_corrs))
clone_rt_corrs <- clone_rt_corrs[, -c('dataset', 'type', 'signature', 'V1')]
clone_cn_dists <- clone_cn_dists[, -c('dataset', 'type', 'signature', 'V1')]
print(colnames(clone_rt_corrs))

# create a mapping of each dataset to a color
# use a for loop to iterate over each dataset
dataset_colors <- c(
    'SA1091' = '#1f77b4',
    'SA1047' = '#ff7f0e',
    'SA1093' = '#2ca02c',
    'SA1049' = '#d62728',
    'SA1096' = '#4C7534',
    'SA1162' = '#661089',
    'SA1050' = '#80BDB6',
    'SA1051' = '#A99B66',
    'SA1052' = '#A43843',
    'SA1053' = '#1D306E',
    'SA1181' = '#B245B2',
    'SA1182' = '#DFE583',
    'SA1184' = '#4F3411',
    'SA530' = '#0F5F5E',
    'SA604' = '#950F24',
    'SA609' = '#3A583C',
    'SA501' = '#757F8D',
    'SA535' = '#EF8151'
)

# create a mapping of each cancer type to a color
type_colors <- c(
    "HGSOC" = "#1D8322",
    "TNBC" = "#AE52D3"
)

# create a mapping of signatures to colors
signature_colors <- c(
    'HRD-Dup' = '#1f77b4',
    'HRD-Del' = '#ff7f0e',
    'FBI' = '#2ca02c',
    'TD' = '#d62728'
)

# plot a complex heatmap of the correlation and distance matrices
# annotate each row with categorial data (dataset, cancer type, signature)
row_ha = rowAnnotation(dataset = datasets_df$dataset, type = datasets_df$type, signature = datasets_df$signature,
                       col = list(dataset = dataset_colors, type = type_colors, signature = signature_colors))
h1 <- Heatmap(
    sample_rt_corrs, name = "Sample RT corr", left_annotation = row_ha, column_title = "Sample RT correlation",
    show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 2)
)
print('Made h1')
h2 <- Heatmap(
    sample_cn_dists, name = "Sample CN dist", left_annotation = row_ha, column_title = "Sample CN distance",
    show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 2)
)
print('Made h2')
# plot complex heatmaps for the clone correlation and distance matrices
clone_row_ha = rowAnnotation(dataset = clone_df$dataset, type = clone_df$type, signature = clone_df$signature,
                             col = list(dataset = dataset_colors, type = type_colors, signature = signature_colors))
h3 <- Heatmap(
    clone_rt_corrs, name = "Clone RT corr", left_annotation = clone_row_ha, column_title = "Clone RT correlation",
    show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.2)
)
print('Made h3')
h4 <- Heatmap(
    clone_cn_dists, name = "Clone CN dist", left_annotation = clone_row_ha, column_title = "Clone CN distance",
    show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = FALSE,
    rect_gp = gpar(col = "white", lwd = 0.2)
)
print('Made h4')

# plot the sample heatmaps
png(args$sample_heatmap, width = 8, height = 4, units = "in", res = 300)
print(ComplexHeatmap::draw(h1 + h2,
                            ht_gap = unit(0.2, "cm"),
                            annotation_legend_side = "bottom",
                            show_heatmap_legend = TRUE))
dev.off()

# plot the clone heatmaps
png(args$clone_heatmap, width = 8, height = 4, units = "in", res = 300)
print(ComplexHeatmap::draw(h3 + h4,
                            ht_gap = unit(0.2, "cm"),
                            annotation_legend_side = "bottom",
                            show_heatmap_legend = TRUE))
dev.off()
