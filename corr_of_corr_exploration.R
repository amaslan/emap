# Generate plots of correlation of correlations for:
#    1. protein vs. mutant - 1 per cluster
#    2. cluster vs. mutant - 1 per protein
#    3. protein-cluster vs. mutant 

library(tidyverse)
library(gplots)
library(RColorBrewer) 
library(rafalib)

######################################
# SPECIFY FILE PATHS, ETC.

# file paths
corr_of_corr = "/Users/annie/emap/corr_of_corr.RData"
output_path = "/Users/annie/emap/20180111/"

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
ppi_only <- filter(sep_no_na, protein %in% ppi)
mut_ppi_only <- filter(ppi_only, mutant %in% mut_data)

# remove grouped clusters
final <- filter(ppi_only, grepl("_GO_", cluster))
#final <- filter(ppi_only, !grepl("_GO_", cluster))

######################################
# PLOT

plot_heatmap <- function(df, id) {
  # create matrix for heatmap plot
  mat <- spread(df, mutant, value)
  rownames(mat) <- mat[,c(1)]
  mat <- mat[,-1]
  mat <- data.matrix(mat)
  
  hmcol <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
  
  par(mar=c(7,4,4,2)+0.1) 
  png(paste(output_path, method, "_", id, ".png", sep=""), width=1600, height=1600)
  gplots::heatmap.2(mat, 
                    trace="none",
                    col=hmcol,
                    margins=c(16,14),
                    key.title=NA,
                    main=id)
  dev.off()
}

# protein vs. mutant - 1 per cluster
for (c in unique(final$cluster_number)) {
    s <- final[which(final$cluster_number == c),]
    s <- s[,c('mutant', 'protein', 'value')]
    plot_heatmap(s, c)
  }

# cluster vs. mutant - 1 per protein
for (p in unique(final$protein)) {
  s <- final[which(final$protein == p),]
  s <- s[,c('mutant', 'cluster', 'value')]
  plot_heatmap(s, p)
}

# protein-cluster vs. mutant
final$pc <- paste(final$protein, "-", final$cluster, sep="")
s <- final[,c('mutant', 'pc', 'value')]
plot_heatmap(s, "all")