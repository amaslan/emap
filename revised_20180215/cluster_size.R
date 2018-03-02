# determine number of genes per cluster for real and random clusters

library(tidyverse)

clusters_real = "/Users/annie/emap/revised_20180215/GO_SLIMS_CLUSTERS_2018-01-17_GO_slims_pearson_complete_clusters.txt"
clusters_rand = "/Users/annie/emap/revised_20180215/RANDOM_CLUSTERS_simulated_cluster_3.txt"

cl_real <- read.table(clusters_real, 
                        header = TRUE,
                        sep="\t",
                        fill=TRUE,
                        quote="")
cl_rand <- read.table(clusters_rand, 
                      header = TRUE,
                      sep="\t",
                      fill=TRUE,
                      quote="")

cl_real_sum <- cl_real %>%
  group_by(cluster) %>%
  dplyr::summarize(num = length(unique(gene_name)))


cl_rand_sum <- cl_rand %>%
  group_by(cluster) %>%
  dplyr::summarize(num = length(unique(gene_name)))

# filter so that the count is for genes that are in the library we screen against

raw_emap = "/Users/annie/emap/June2016_Gsp1_E-MAP_data.RData"
load(raw_emap)
all <- e.map
emap_genes = unique(all$library_gene_name)

cl_real_sum_emap <- cl_real %>%
  filter(gene_name %in% emap_genes) %>%
  group_by(cluster) %>%
  dplyr::summarize(num = length(unique(gene_name)))

cl_rand_sum_emap <- cl_rand %>%
  filter(gene_name %in% emap_genes) %>%
  group_by(cluster) %>%
  dplyr::summarize(num = length(unique(gene_name)))

