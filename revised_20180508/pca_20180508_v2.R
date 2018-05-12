# principal component analysis
# alter first section to specify: (1) output path, (2) corr of corr file,
# (3) subset of mutations, (4) subset of genes
# known bug: peroxisome PCA with all genes has SVD error

library(tidyverse)
library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)

#################################################################
### Specify file paths for output and which correlation of correlations file
output_path = "/Users/annie/emap/revised_20180508/pca_results/"
corr_of_corr = "/Users/annie/emap/revised_20180508/20180507_corr_of_corr_filter_FDR.RData"

load(corr_of_corr)


#### only keep cluster + geneB combinations where at least one of the mutants has a very high correlation
### i.e. absolute value larger than 0.7
### but don't keep if the abs(corr) in WT negative control (GSP1-NAT) is larger than 0.3
gene_cluster_pairs_with_strong_WT_corr <- combined_filtered_data %>% 
  filter(geneA == "GSP1-NAT" & abs(corr) > 0.3) %>% 
  mutate("uniq_gene_name" = str_c(geneB, gene_name, sep = " ")) %>% 
  select(geneA, uniq_gene_name, cluster)
strong_corr_subset <- combined_filtered_data %>% 
  group_by(geneB, cluster) %>%
  mutate("per_mutant_max_corr" = max(corr)) %>% 
  filter(per_mutant_max_corr > 0.7) %>% 
  ungroup() %>% 
  mutate("uniq_gene_name" = str_c(geneB, gene_name, sep = " "))

strong_corr_subset <- anti_join(strong_corr_subset, 
                                gene_cluster_pairs_with_strong_WT_corr, by = c("uniq_gene_name", "cluster"))


pca_input <- strong_corr_subset %>% spread(geneA, corr)

#################################################################
# PCA by gene
# can we distinguish the different sub GO categories from
# hierarchical clustering based on PCA?

# look for sub GO categories

cls <- c(
  "cell cycle",
  "chromatin",
  "cellular response to DNA damage stimulus",
  "cytoskeleton",
  "DNA",
  "endomembrane",
  "endoplasmic",
  "ER",
  "histone",
  "lipid",
  "mitochondri",
  "organelle",
  "protein",
  "transcription",
  "whole_library"
)

for (i in 1:length(cls)) { 
  pca_input_cluster <- pca_input %>% filter(grepl(cls[i], cluster))
  
  # compute principal components
  to_plot <- pca_input_cluster[,-c(1:11)]
  to_plot[is.na(to_plot)] <- 0
  gene_pca <- prcomp(to_plot, center=TRUE, scale=TRUE)

  # print(gene_pca)
  # summary(gene_pca)
  
  # plot of variances (y-axis) associated with the PCs (x-axis)
  # plot(gene_pca, type = "l")
  
  gr <- as.factor(pca_input_cluster$cluster)
  
  # project data on the first two PCs
  g <- ggbiplot(gene_pca, obs.scale = 1, var.scale = 1, 
                groups = gr, ellipse = TRUE) + 
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal', 
          legend.position = 'top',
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) 
  
  ggsave(filename=paste(output_path, "all_", cls[i], "_pca", ".png", sep=""), width = 10, height = 10)
  }

# #################################################################
# # repeat PCA analysis by mutant rather than by gene
# # can we distinguish the different GO categories based on PCA?
# 
# pca_input <- combined_filtered_data_nona %>% spread(gene_name, corr)
# 
# 
# # long <- combined_filtered_data_nona %>% 
# #   tibble::rowid_to_column() %>% 
# #   gather(gene_name, corr, -rowid)
# # 
# # pca_input <- long %>% 
# #   spread(gene_name, corr)
# 
# #pca_input <- combined_filtered_data_nona %>% spread(gene_name, corr)
# 
# for (i in 1:length(cls)) {
#   pca_input_cluster <- pca_input %>% filter(grepl(cls[i], cluster))
# 
#   # compute principal components
#   gene_pca <- prcomp(pca_input_cluster[,-c(1,2)],center=TRUE, scale=TRUE)
#   #print(gene_pca)
#   #summary(gene_pca)
# 
#   # plot of variances (y-axis) associated with the PCs (x-axis)
#   plot(gene_pca, type = "l")
# 
#   gr <- as.factor(pca_input_cluster$cluster)
# 
#   # project data on the first two PCs
#   g <- ggbiplot(gene_pca, obs.scale = 1, var.scale = 1,
#                 groups = gr, ellipse = TRUE) +
#     scale_color_discrete(name = '') +
#     theme(legend.direction = 'horizontal',
#           legend.position = 'top',
#           panel.background = element_blank(),
#           axis.line = element_line(colour = "black"))
# 
#   ggsave(filename=paste(output_path, "all_mutant_", cls[i], "_pca", ".png", sep=""), width = 10, height = 10)
# }
# 
# #################################################################
# # by gene
# # can we distinguish the different GO categories based on PCA for most important categories:
# # used old data that has merged info for these!!!
# 
# pca_input_imp <- pca_input %>% filter((cluster == 'chromatin') 
#                                       | (cluster == 'cell cycle') 
#                                       | (cluster == 'transcription and mRNA processing') 
#                                       | (cluster == 'Golgi and ER'))
#   
# gene_pca <- prcomp(pca_input_imp[,-c(1,2)], center=TRUE, scale=TRUE)
#   
# gr <- as.factor(pca_input_imp$cluster)
#   
# # project data on the first two PCs
# g <- ggbiplot(gene_pca, obs.scale = 1, var.scale = 1, 
#                 groups = gr, ellipse = TRUE) + 
#     scale_color_discrete(name = '') +
#     theme(legend.direction = 'horizontal', 
#           legend.position = 'top',
#           panel.background = element_blank(), 
#           axis.line = element_line(colour = "black")) #+
#   
# ggsave(filename=paste(output_path, "all_imp_pca", ".png", sep=""), width = 10, height = 10)
# 


# ### OLD
# combined_filtered_data <- combined_filtered_data[,c("corr", "geneA", "gene_name", "cluster")]
# 
# # convert value to numeric and remove na
# combined_filtered_data$corr <- as.numeric(as.character(combined_filtered_data$corr))
# combined_filtered_data_nona <- filter(combined_filtered_data, !is.na(corr))
# 
# # additional filter to include strong mutants only v2 with more included
# strong_muts1 <- c("D79A", "D79S", "H141E", "H141R", "K101R", "R108L", "R108Q", "R108Y",
#                   "R112A", "R112S", "R78K", "T34E", "T34G", "T34Q") 
# strong_muts2 <- c("D79A", "D79S", "H141E", "H141R", "K101R", "R108L", "R108Q", "R108Y", 
#                   "R112A", "R112S", "R78K", "T34E", "T34G", "T34Q", "R108G", "H141I",
#                   "Y148I", "R108S", "R108A", "T34N","Y157A")
# combined_filtered_data_nona <- combined_filtered_data_nona[which((combined_filtered_data_nona$geneA %in% strong_muts2)),]
# 
# # additional filter to only include genes for which we have data
# #ppi = c("MOG1", "KAP104", "MTR10", "MSN5", "LOS1", "CRM1", "PSE1", "KAP95", "SRP1", "YRB1", "YRB2", "SRM1", "RNA1", "NTF2", "CSE1", "NUP60")
# 
# 
# #combined_filtered_data_nona <- combined_filtered_data_nona[which((combined_filtered_data_nona$gene_name %in% ppi)),]
# 
# #dup_genes = c("ACT1", "CDC4", "DED81", "DSN1", "GIM4", "LSM3", "MEX67", "MPS1", "NUP57", "PRP4", "RVS167", "SEC22", "SMC1", "SMC3", "STU1", "WBP1", "YHC1", "YNL181W")
# #combined_filtered_data_nona <- combined_filtered_data_nona[which(!(combined_filtered_data_nona$gene_name %in% dup_genes)),]
# 
# # i=1
# # for (i in 1:dim(combined_filtered_data_nona)[1]) {
# #   combined_filtered_data_nona[i,]$gene_name <- paste(combined_filtered_data_nona[i,]$gene_name, "-", i)
# # }
# 
# combined_filtered_data_nona$ID <- seq.int(nrow(combined_filtered_data_nona))
# 
# pca_input <- combined_filtered_data_nona %>% spread(geneA, corr)
# pca_input <- pca_input[,-c(3)]
# 
# #long <- combined_filtered_data_nona %>% 
# #tibble::rowid_to_column() %>% 
# #gather(geneA, corr, -rowid)
# 
# #pca_input <- long %>% 
# #spread(geneA, corr)