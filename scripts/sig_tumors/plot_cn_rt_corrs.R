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
parser$add_argument("--heatmap", default=NULL, type="character",
                    help="heatmap plot")
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

print(sample_rt_corrs$dataset)
print(sample_cn_dists$dataset)

# check that the order of the samples in the correlation and distance matrices match
stopifnot(all(sample_rt_corrs$dataset == sample_cn_dists$dataset))

# sort the rows of datasets_df to match the order of the rows in the correlation and distance matrices
print(datasets_df)
datasets_df <- datasets_df[match(sample_rt_corrs$dataset, datasets_df$dataset),]
print(datasets_df)

# remove the dataset and V1 columns from the correlation and distance matrices
sample_rt_corrs <- sample_rt_corrs[, -c('dataset', 'V1')]
sample_cn_dists <- sample_cn_dists[, -c('dataset', 'V1')]

# plot a complex heatmap of the correlation and distance matrices
# annotate each row with categorial data (dataset, cancer type, signature)
row_ha = rowAnnotation(dataset = datasets_df$dataset, type = datasets_df$type, signature = datasets_df$signature)
h1 <- Heatmap(sample_rt_corrs, name = "RT corr", left_annotation = row_ha)
print('Made h1')
h2 <- Heatmap(sample_cn_dists, name = "CN dist", left_annotation = row_ha)
print('Made h2')

png(args$heatmap, width = 25, height = 12, units = "in", res = 300)
print(ComplexHeatmap::draw(h1 + h2,
                            ht_gap = unit(0.6, "cm"),
                            # heatmap_legend_side = "bottom",
                            annotation_legend_side = "bottom",
                            show_heatmap_legend = TRUE))
dev.off()
