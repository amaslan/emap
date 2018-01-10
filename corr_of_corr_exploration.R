library(ggplot2)
library(dplyr)
library(tidyr)

# file paths
corr_of_corr = "/Users/annie/emap/corr_of_corr.RData"
output_path = "/Users/annie/emap/20180109/"

# mutants for which we have data
mut_data <- c("GSP1-NAT", "T34E", "R108L", "H141V", "Q147E")

# PPI partners 
ppi = c("MOG1", "KAP104", "MTR10", "MSN5", "LOS1", "CRM1", "PSE1", "KAP95", "SRP1", "YRB1", "YRB2", "SRM1", "RNA1", "NTF2", "CSE1")

load(corr_of_corr)

######################################
# FILTER

# look only at method corr_of_corr_no_na_no_mut
corr_of_corr_no_na_no_mut <- result$corr_of_corr[which(result$corr_of_corr$method == "corr_of_corr_no_na_no_mut"),]

# separate partner column into protein and ORF
sep <- corr_of_corr_no_na_no_mut %>% 
  separate(partner, into=c("protein", "ORF"), sep=" - ", remove=FALSE)

write.table(sep, file=paste(output_path, "sep.csv"), quote=FALSE, sep=",", row.names=FALSE)

# filter to only include PPI partners and mutants for which we have data
ppi_only <- filter(sep, protein %in% ppi)
mut_ppi_only <- filter(ppi_only, mutant %in% mut_data)

# convert value to numeric and remove na
mut_ppi_only$value <- as.numeric(mut_ppi_only$value)
mut_ppi_no_na_only <- filter(mut_ppi_only, !is.na(value))

# remove grouped clusters
final <- filter(mut_ppi_no_na_only, grepl("_GO_", cluster))


######################################
# PLOT

# protein vs. mutant - 1 per cluster
for (c in unique(final$cluster_number)) {
  ggplot(data = final[which(final$cluster_number == c),], aes(x=mutant, y=protein, fill=value)) + 
        geom_tile() +
        labs(title=final[which(final$cluster_number == c),]$cluster)
  ggsave(filename=paste(output_path, c, "_protein_v_mutant", ".png", sep=""), width = 10, height = 10)
  }

# cluster vs. mutant - 1 per protin
for (p in unique(final$protein)) {
  ggplot(data = final[which(final$protein == p),], aes(x=mutant, y=cluster, fill=value)) + 
    geom_tile() +
    labs(title=p)
  ggsave(filename=paste(output_path, p, "_cluster_v_mutant", ".png", sep=""), width = 10, height = 10)
}

write.table(final, file=paste(output_path, "final.csv"), quote=FALSE, sep=",", row.names=FALSE)


