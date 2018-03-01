# Generate plots of correlation of correlations for:
#    1. protein vs. mutant - 1 per cluster
#    2. cluster vs. mutant - 1 per protein
#    3. protein-cluster vs. mutant 

library(tidyverse)
library(gplots)
library(RColorBrewer) 
library(rafalib)
library(reshape2)

######################################
# SPECIFY FILE PATHS, ETC.

# file paths
#corr_of_corr = "/Users/annie/emap/revised_20180215/random_corr_of_corr.RData"
corr_of_corr = "/Users/annie/emap/revised_20180215/GO_slims_corr_of_corr.RData"
output_path = "/Users/annie/emap/revised_20180215/20180228/"

# PPI partners - added NUP60
ppi = c("MOG1", "KAP104", "MTR10", "MSN5", "LOS1", "CRM1", "PSE1", "KAP95", "SRP1", "YRB1", "YRB2", "SRM1", "RNA1", "NTF2", "CSE1", "NUP60")

# method
# look only at method corr_of_corr_no_na_no_mut
# will also be interested in dot_product_no_na_no_mut
method = "corr_of_corr_no_mut"

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

# remove grouped clusters
# final <- filter(ppi_only, grepl("_GO_", cluster))
final <- sep_no_na

# just look at strong muts
#strong_muts1 <- c("D79A", "D79S", "H141E", "H141R", "K101R", "R108L", "R108Q", "R108Y",
                 # "R112A", "R112S", "R78K", "T34E", "T34G", "T34Q") 
mut_filter <- c("R108L", "R108Q", "R108Y",
                 "T34E", "T34G", "T34Q") 

#strong_muts1 <- c("D79A", "D79S", "H141E", "H141R",
                  #"R112A", "R112S") 

# mut_filter <- c("R108L", "R108Q", "R108Y", 
#                 "R108G",
#                 "R108S", "R108A")

# mut_filter <- c("D79A", "D79S", "H141E", "H141R", "K101R", "R108L", "R108Q", "R108Y", 
#                   "R112A", "R112S", "R78K", "T34E", "T34G", "T34Q", "R108G", "H141I",
#                   "Y148I", "R108S", "R108A", "T34N","Y157A")

final <- filter(final, mutant %in% mut_filter)
#final <- filter(final, grepl('T34', mutant))

# partner vs. mutant - 1 per cluster
i=1
m <- c()
for (c in unique(final$cluster)) {
    s <- final[which(final$cluster == c),]
    s <- s[,c('mutant', 'partner', 'value')]
    mat <- spread(s, partner, value)
    rownames(mat) <- mat[,c(1)]
    mat <- mat[,-1]
    mat <- data.matrix(mat)
    d <- dist(mat)
    tmp <- data.frame(t(combn(rownames(mat),2)), as.numeric(d), c)
    m <- rbind(m, tmp)
    i = i+1
}

names(m) <- c("p1", "p2", "distance", "cluster")

m$pair <- paste(m$p1, "-", m$p2)

mat <- spread(m, cluster, distance)
rownames(mat) <- mat[,c(3)]
mat <- mat[,!names(mat) %in% c('p1', 'p2', 'pair')]
mat <- data.matrix(mat)
d <- dist(mat)

hc <- hclust(d)
plot(hc)

fit <- cmdscale(d, eig=TRUE, k=2)
fit

x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="coordinate 1", ylab='coordinate 2', main='MDS', type='n')
text(x, y, labels = row.names(mat), cex=.5)


f <- data.frame(t(combn(rownames(mat),2)), as.numeric(d))
names(f) <- c("pair1", "pair2", "distance")


ggplot(f, aes(x=pair1, y=distance)) +
  geom_boxplot(colour = "#1F3552", fill = "#4271AE",
               size = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



dev.off()
hmcol <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
gplots::heatmap.2(as.matrix(d), 
                  trace="none",
                  col=hmcol,
                  #margins=c(1,1),
                  key.title=NA,
                  main=id)
