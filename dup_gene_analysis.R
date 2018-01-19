# analyze correlation of correlations differences across experiments
# for the 18 genes with 2-6 experiments. for subsequent analysis
# we take the correlation of correlations value with the largest
# absolute value

library(tidyverse)

# determine distribution of correlation of correlations
all_rds = "/Users/annie/emap/all.rds"
output_path = "/Users/annie/emap/20180118/"

all = readRDS(all_rds)

# remove duplicates that arise because multiple interface scores of 0
#all <- all[!duplicated(all),]
# remove complex and interface scores info so don't have duplicate correlation of correlation entries
all <- all[,!(names(all) == "Complex")]
all <- all[,!(names(all) == "interface.score")]
all <- all[!duplicated(all),]
# remove NA correlation of correlation values
all$corr.value <- as.numeric(all$corr.value)
all_no_na <- filter(all, !is.na(corr.value))

dups <- aggregate(all_no_na$gene_mut, list(all_no_na$cluster, all_no_na$protein), length)
dups <- dups[which(dups$x > length(unique(all_no_na$mutant))),]
dup_genes <- unique(dups$Group.2)

all_no_na_dups <- all_no_na[(which(all_no_na$protein %in% dup_genes)),]

# remove data for merged clusters
all_no_na_dups_no_merge <- filter(all_no_na_dups, grepl("_GO_", cluster))

by_gene_mut_cluster <- group_by(all_no_na_dups_no_merge, protein, mutant, cluster)

by_gene_mut_cluster_dup_summary <- dplyr::summarize(by_gene_mut_cluster, 
          total.count = n(), 
          min.val = min(corr.value), 
          max.val = max(corr.value),
          diff = max(corr.value) - min(corr.value),
          mean.val = mean(corr.value))

final <- filter(by_gene_mut_cluster_dup_summary, total.count > 1)

by_gene <- group_by(final, protein)
summary_gene <- dplyr::summarize(by_gene, 
                              mean.diff = mean(diff),
                              median.diff = median(diff),
                              min.diff = min(diff),
                              max.diff = max(diff))

write.table(summary_gene, file=paste(output_path, "duplicate_summary_by_gene.csv"), quote=FALSE, sep=",", row.names=FALSE)

by_mutant <- group_by(final, mutant)
summary_mutant <- dplyr::summarize(by_mutant, 
                                 mean.diff = mean(diff),
                                 median.diff = median(diff),
                                 min.diff = min(diff),
                                 max.diff = max(diff))

by_cluster <- group_by(final, cluster)
summary_cluster <- dplyr::summarize(by_cluster, 
                                   mean.diff = mean(diff),
                                   median.diff = median(diff),
                                   min.diff = min(diff),
                                   max.diff = max(diff))

ggplot(data=final, aes(x=diff)) +
  geom_histogram(aes(y=..density..), color="black", fill="white", bins = 100) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(final$diff, na.rm=TRUE)), color="blue", linetype="dashed", size=1) +
  labs(x="max - min corr of corr") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data=final, aes(x=mean.val, y=diff)) +
  geom_point(size=0.001) +
  labs(x="mean corr of corr", y="max - min corr of corr") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data=final, aes(x=max.val, y=diff)) +
  geom_point(size=0.001) +
  labs(x="max corr of corr", y="max - min corr of corr") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data=final, aes(x=min.val, y=diff)) +
  geom_point(size=0.001) +
  labs(x="min corr of corr", y="max - min corr of corr") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data=final, aes(x=min.val, y=max.val)) +
  geom_point(size=0.001) +
  labs(x="min corr of corr", y="max corr of corr") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

mean(final$diff)
quantile(final$diff)
