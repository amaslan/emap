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
output_path = "/Users/annie/emap/revised_20180215/20180227/"
corr_of_corr = "/Users/annie/emap/revised_20180215/GO_slims_corr_of_corr.RData"
# corr_of_corr = "/Users/annie/emap/revised_20180215/random_corr_of_corr.RData"
# corr_of_corr = "/Users/annie/emap/corr_of_corr.RData"

# method for correlation of correlations file - depends on if old or new corr of corr
method = "corr_of_corr_no_mut"
# method = "corr_of_corr_no_na_no_mut"

# result$corr_of_corr: (1) method (2) value (3) mutant (4) partner (5) cluster_number (6) cluster
# e.g.: 1      corr_of_corr_all  -0.423138711324343  K129F CIN1 - YOR349W             10 Golgi and ER_GO_2
load(corr_of_corr)
result_m <- result$corr_of_corr[which(result$corr_of_corr$method == method),]
sep <- result_m %>% 
  separate(partner, into=c("protein", "ORF"), sep=" - ", remove=FALSE)
sep_filtered <- sep[,c("value", "mutant", "partner", "cluster")]
sep_filtered$value <- as.numeric(as.character(sep_filtered$value))
sep_filtered_no_na <- filter(sep_filtered, !is.na(value))

# additional filter to include strong mutants only v2 with more included
strong_muts1 <- c("D79A", "D79S", "H141E", "H141R", "K101R", "R108L", "R108Q", "R108Y",
                  "R112A", "R112S", "R78K", "T34E", "T34G", "T34Q") 
strong_muts2 <- c("D79A", "D79S", "H141E", "H141R", "K101R", "R108L", "R108Q", "R108Y", 
                  "R112A", "R112S", "R78K", "T34E", "T34G", "T34Q", "R108G", "H141I",
                  "Y148I", "R108S", "R108A", "T34N","Y157A")
sep_filtered_no_na <- sep_filtered_no_na[which((sep_filtered_no_na$mutant %in% strong_muts2)),]

# additional filter to only include genes for which we have data
partners_with_data <- c("RNA1 - YMR235C_TSQ172", "SRM1 - YGL097W_TSQ958", 
                        "MOG1 - YJR074W", "YRB1 - YDR002W_TSQ582", "PSE1 - YMR308C_TSQ683", "NUP60 - YAR002W")

sep_filtered_no_na <- sep_filtered_no_na[which((sep_filtered_no_na$partner %in% partners_with_data)),]

pca_input <- sep_filtered_no_na %>% spread(mutant, value)

#################################################################
# PCA by gene
# can we distinguish the different sub GO categories from
# hierarchical clustering based on PCA?

# removed merged
cls <- c("budding", "cell cycle", "chromatin", 
         "cytoskeleton", "Golgi", "lipids", 
         "metabolic", "mitochondrion", 
         "nuclear", "peroxisome", "ribosome",
         "transcription", "vacuole")

for (i in 1:length(cls)) { 
  pca_input_cluster <- pca_input %>% filter(grepl(cls[i], cluster))
  
  # compute principal components
  gene_pca <- prcomp(pca_input_cluster[,-c(1,2)], center=TRUE, scale=TRUE)
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
          axis.line = element_line(colour = "black")) #+
  
  ggsave(filename=paste(output_path, "all_", cls[i], "_pca", ".png", sep=""), width = 10, height = 10)
  }

#################################################################
# repeat PCA analysis by mutant rather than by gene
# can we distinguish the different GO categories based on PCA?

pca_input <- sep_filtered_no_na %>% spread(partner, value)

for (i in 1:length(cls)) {
  pca_input_cluster <- pca_input %>% filter(grepl(cls[i], cluster))

  # compute principal components
  gene_pca <- prcomp(pca_input_cluster[,-c(1,2)],center=TRUE, scale=TRUE)
  #print(gene_pca)
  #summary(gene_pca)

  # plot of variances (y-axis) associated with the PCs (x-axis)
  plot(gene_pca, type = "l")

  gr <- as.factor(pca_input_cluster$cluster)

  # project data on the first two PCs
  g <- ggbiplot(gene_pca, obs.scale = 1, var.scale = 1,
                groups = gr, ellipse = TRUE) +
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal',
          legend.position = 'top',
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))

  ggsave(filename=paste(output_path, "all_mutant_", cls[i], "_pca", ".png", sep=""), width = 10, height = 10)
}

#################################################################
# by gene
# can we distinguish the different GO categories based on PCA for most important categories:
# used old data that has merged info for these!!!

pca_input_imp <- pca_input %>% filter((cluster == 'chromatin') 
                                      | (cluster == 'cell cycle') 
                                      | (cluster == 'transcription and mRNA processing') 
                                      | (cluster == 'Golgi and ER'))
  
gene_pca <- prcomp(pca_input_imp[,-c(1,2)], center=TRUE, scale=TRUE)
  
gr <- as.factor(pca_input_imp$cluster)
  
# project data on the first two PCs
g <- ggbiplot(gene_pca, obs.scale = 1, var.scale = 1, 
                groups = gr, ellipse = TRUE) + 
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal', 
          legend.position = 'top',
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) #+
  
ggsave(filename=paste(output_path, "all_imp_pca", ".png", sep=""), width = 10, height = 10)