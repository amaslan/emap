# Generate plots of correlation of correlations for:
#    1. protein vs. mutant - 1 per cluster
#    2. cluster vs. mutant - 1 per protein
#    3. protein-cluster vs. mutant 

library(ggplot2)
library(dplyr)
library(tidyr)

######################################
# SPECIFY FILE PATHS, ETC.

# file paths
corr_of_corr = "/Users/annie/emap/corr_of_corr.RData"
output_path = "/Users/annie/emap/20180110/"

# mutants for which we have data
mut_data <- c("GSP1-NAT", "T34E", "R108L", "H141V", "Q147E")

# PPI partners 
ppi = c("MOG1", "KAP104", "MTR10", "MSN5", "LOS1", "CRM1", "PSE1", "KAP95", "SRP1", "YRB1", "YRB2", "SRM1", "RNA1", "NTF2", "CSE1")

# method
# look only at method corr_of_corr_no_na_no_mut
# will also be interested in dot_product_no_na_no_mut
method = "corr_of_corr_no_na_no_mut"

######################################
# LOAD, CLEAN, AND FILTER

load(corr_of_corr)

result_m <- result$corr_of_corr[which(result$corr_of_corr$method == method),]

# separate partner column into protein and ORF
sep <- result_m %>% 
  separate(partner, into=c("protein", "ORF"), sep=" - ", remove=FALSE)

# convert value to numeric and remove na
sep$value <- as.numeric(sep$value)
sep_no_na <- filter(sep, !is.na(value))

# save data as csv for input to python code for MDS
write.table(sep_no_na, file=paste(output_path, method, "_sep_no_na.csv"), quote=FALSE, sep=",", row.names=FALSE)

# filter to only include PPI partners and mutants for which we have data
ppi_only <- filter(sep, protein %in% ppi)
mut_ppi_only <- filter(ppi_only, mutant %in% mut_data)

# remove grouped clusters
# final <- filter(mut_ppi_only, grepl("_GO_", cluster))
final <- filter(ppi_only, grepl("_GO_", cluster))



######################################
# PLOT

# protein vs. mutant - 1 per cluster
for (c in unique(final$cluster_number)) {
  ggplot(data = final[which(final$cluster_number == c),], aes(x=mutant, y=protein, fill=value)) + 
        geom_tile() +
        labs(title=final[which(final$cluster_number == c),]$cluster) +
        scale_fill_gradient(low = "purple", high = "pink") +
        theme(axis.text.x = element_text(angle = 90, hjust=1)) +
        coord_equal()
  ggsave(filename=paste(output_path, c, "_", method, "_protein_v_mutant", ".png", sep=""), width = 10, height = 10)
  }

# cluster vs. mutant - 1 per protein
for (p in unique(final$protein)) {
  ggplot(data = final[which(final$protein == p),], aes(x=mutant, y=cluster, fill=value)) + 
    geom_tile() +
    labs(title=p) +
    scale_fill_gradient(low = "purple", high = "pink") +
    theme(axis.text.x = element_text(angle = 90, hjust=1)) +
    coord_equal()
  ggsave(filename=paste(output_path, p, "_", method, "_cluster_v_mutant", ".png", sep=""), width = 10, height = 10)
}

# protein-cluster vs. mutant
final$pc <- paste(final$protein, "-", final$cluster, sep="")
ggplot(data = final, aes(x=mutant, y=pc, fill=value)) + 
  geom_tile() +
  labs(title=all) +
  scale_fill_gradient(low = "purple", high = "pink") +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  coord_equal()
ggsave(filename=paste(output_path, method, "_protein-cluster_v_mutant", ".png", sep=""), width = 20, height = 20)



