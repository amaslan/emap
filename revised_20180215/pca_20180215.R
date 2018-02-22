# principal component analysis

library(tidyverse)
library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)


output_path = "/Users/annie/emap/revised_20180215/20180215/"

corr_of_corr = "/Users/annie/emap/revised_20180215/GO_slims_corr_of_corr.RData"

# method for correlation of correlations file
method = "corr_of_corr_no_mut"

# result$corr_of_corr: (1) method (2) value (3) mutant (4) partner (5) cluster_number (6) cluster
# e.g.: 1      corr_of_corr_all  -0.423138711324343  K129F CIN1 - YOR349W             10 Golgi and ER_GO_2
load(corr_of_corr)
result_m <- result$corr_of_corr[which(result$corr_of_corr$method == method),]
sep <- result_m %>% 
  separate(partner, into=c("protein", "ORF"), sep=" - ", remove=FALSE)


sep_filtered <- sep[,c("value", "mutant", "protein", "cluster")]
sep_filtered$value <- as.numeric(as.character(sep_filtered$value))
sep_filtered_no_na <- filter(sep_filtered, !is.na(value))

sep_filtered_no_na <- sep_filtered_no_na[which(!(sep_filtered_no_na$protein %in% dup_genes)),]

# additional filter to include strong mutants only v2 with more included
sep_filtered_no_na <- sep_filtered_no_na[which((sep_filtered_no_na$mutant %in% strong_muts2)),]

# additional filter to only include genes for which we have data
genes_with_data <- c("RNA1", "SRM1", "MOG1", "YRB1", "PSE1")
#sep_filtered_no_na <- sep_filtered_no_na[which((sep_filtered_no_na$protein %in% genes_with_data)),]


pca_input <- sep_filtered_no_na %>% spread(mutant, value)

#################################################################
# # can we distinguish the complex membership based on PCA?
# # check for each cluster
# 
# complex_def = "/Users/annie/emap/CYC2008_complex.tab.txt"
# complex <- read.table(complex_def,
#                       header = TRUE,
#                       sep="\t",
#                       col.names=c("ORF", "Name", "Complex", "PubMed_id", "Method", "GO_id",
#                                   "GO_term", "Jaccard_Index"),
#                       fill=TRUE,
#                       quote="")
# 
# # just care about ORF, Name, Complex from complex table
# keep = c("ORF", "Name", "Complex")
# complex_short <- complex[names(complex) %in% keep]
# 
# pca_input_2 <- merge(complex_short, pca_input, by.x="Name", by.y="protein")
# 
# for (c in unique(pca_input_2$cluster)) {
#   pca_input_cluster <- pca_input_2 %>% filter(cluster == c)
#   
#   # compute principal components
#   gene_pca <- prcomp(pca_input_cluster[,-c(1,2,3,4)], center=TRUE, scale=TRUE)
#   #print(gene_pca)
#   #summary(gene_pca)
#   
#   # plot of variances (y-axis) associated with the PCs (x-axis)
#   #plot(gene_pca, type = "l")
#   
#   gr <- as.factor(pca_input_cluster$Complex)
#   
#   # project data on the first two PCs
#   g <- ggbiplot(gene_pca, obs.scale = 1, var.scale = 1, 
#                 groups = gr, ellipse = TRUE) + 
#     #circle = TRUE) +
#     scale_color_discrete(name = '') +
#     theme(legend.direction = 'horizontal', 
#           legend.position = 'none',
#           panel.background = element_blank(), 
#           axis.line = element_line(colour = "black"))
#   #print(g)
#   
#   ggsave(filename=paste(output_path, c, "_pca", ".png", sep=""), width = 10, height = 10)
# }
# 

#################################################################
# by gene
# can we distinguish the different GO categories based on PCA?

#cls <- unique(pca_input$cluster)
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
  #print(gene_pca)
  #summary(gene_pca)
  
  # plot of variances (y-axis) associated with the PCs (x-axis)
  #plot(gene_pca, type = "l")
  
  gr <- as.factor(pca_input_cluster$cluster)
  
  # project data on the first two PCs
  g <- ggbiplot(gene_pca, obs.scale = 1, var.scale = 1, 
                groups = gr, ellipse = TRUE) + 
    #circle = TRUE) +
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal', 
          legend.position = 'top',
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black")) #+
    #ylim(-3, 3)
  #print(g)
  
  ggsave(filename=paste(output_path, "all_", cls[i], "_pca", ".png", sep=""), width = 10, height = 10)
  
  }

#################################################################
# repeat analysis by mutant rather than by gene
# can we distinguish the different GO categories based on PCA?
# just did with merged

pca_input <- sep_filtered_no_na %>% spread(protein, value)

cls <- c("budding", "cell cycle", "chromatin",
         "cytoskeleton", "Golgi", "lipids",
         "merged", "metabolic", "mitochondrion",
         "nuclear", "peroxisome", "ribosome",
         "transcription", "vacuole")

for (i in 1:length(cls)) {
  pca_input_cluster <- pca_input %>% filter(grepl(cls[i], cluster))
  #pca_input_cluster <- pca_input_cluster[colSums(!is.na(pca_input_cluster)) > 0]
  #pca_input_cluster <- na.omit(pca_input_cluster)
    #na.omit(pca_input_cluster)

  # compute principal components
  gene_pca <- prcomp(pca_input_cluster[,-c(1,2)],center=TRUE, scale=TRUE)
  print(gene_pca)
  summary(gene_pca)

  # plot of variances (y-axis) associated with the PCs (x-axis)
  plot(gene_pca, type = "l")

  gr <- as.factor(pca_input_cluster$cluster)

  # project data on the first two PCs
  g <- ggbiplot(gene_pca, obs.scale = 1, var.scale = 1,
                groups = gr, ellipse = TRUE) +
    #circle = TRUE) +
    scale_color_discrete(name = '') +
    theme(legend.direction = 'horizontal',
          legend.position = 'top',
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
  #print(g)

  ggsave(filename=paste(output_path, "all_mutant_", cls[i], "_pca", ".png", sep=""), width = 10, height = 10)

}


