# analyze correlation of correlations differences across experiments
# for the 18 genes with 2-6 experiments. for subsequent analysis
# we take the correlation of correlations value with the largest
# absolute value

library(tidyverse)

# determine distribution of correlation of correlations
all_rds = "/Users/annie/emap/all.rds"
output_path = "/Users/annie/emap/20180118/"

all = readRDS(all_rds)

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

############################################################################
# show what the different experiments are, e.g. DAMP, temperature
# use strongest corr of corr unless check to see if construct with greatest strength
# varies by mutant then keep both and treat as different genes

output_path = "/Users/annie/emap/20180126/"

# start with all_no_na_dups which also includes merged clusters
for (p in dup_genes) {
  all_no_na_p <- all_no_na_dups %>% filter(protein == p)
  ggplot(all_no_na_p, aes(x=partner, y=corr.value)) +
    geom_boxplot(colour = "#1F3552", fill = "#4271AE",
                 size = 1) +
    #geom_jitter() +
    labs(title=c) +
    theme_bw() +
    geom_hline(yintercept = 0, colour="red") +
    theme(panel.grid.major = element_line(colour = "#d3d3d3"),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
            text=element_text(family = "Tahoma"),
            axis.title = element_text(face="bold"),
            axis.text.x = element_text(colour="black", size = 11),
            axis.text.y = element_text(colour="black", size = 9),
            axis.line = element_line(size=0.5, colour = "black"))
    ggsave(filename=paste(output_path, p, "_duplicate_boxplot", ".png", sep=""), width = 10, height = 10)
}


# top$partner has the constructs with the max(mean abs(corr.value)) for each protein
by_construct <- all_no_na_dups %>%
  group_by(protein, partner)
summary_construct <- dplyr::summarize(by_construct, 
                                      med.val = median(abs(corr.value)),
                                      mean.val = mean(abs(corr.value)))
top <- 
  summary_construct %>% 
  group_by(protein) %>% 
  top_n(1, (mean.val))
  #top_n(1, (med.val))

# get what the corr.value is for each gene_mut_cluster for top constructs
all_no_na_dups_top <- filter(all_no_na_dups, partner %in% top$partner)
top_by_gene_mut_cluster <- group_by(all_no_na_dups_top, gene_mut_cluster)
summary_top_gene_mut_cluster <- dplyr::summarize(top_by_gene_mut_cluster, 
                                              top.val = abs(corr.value))

# calculate delta = abs(max partner corr.value) - abs(each other construct corr.value)
all_no_na_dups_final <- merge(all_no_na_dups, summary_top_gene_mut_cluster, by="gene_mut_cluster", all=TRUE)
all_no_na_dups_final$delta <- all_no_na_dups_final$top.val - abs(all_no_na_dups_final$corr.value)

# plot deltas for all
ggplot(data=all_no_na_dups_final, aes(x=delta)) +
  geom_histogram(aes(y=..density..), color="black", fill="white", bins=100) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(all_no_na_dups_final$delta)), color="blue", linetype="dashed", size=1) +
  labs(x="abs(max partner corr.value) - abs(each other construct corr.value)") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

# plot deltas by protein
for (p in unique(all_no_na_dups_final$protein)) {
  sub <- all_no_na_dups_final[which(all_no_na_dups_final$protein == p & !(all_no_na_dups_final$partner %in% top$partner)),]
  ggplot(data=sub, aes(x=delta)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins=100) +
    geom_density(alpha=.2, fill="#FF6666") +
    geom_vline(aes(xintercept=mean(sub$delta, na.rm=TRUE)), color="blue", linetype="dashed", size=1) +
    labs(x="abs(max partner corr.value) - abs(each other construct corr.value)", title=p) +
    theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(filename=paste(output_path, p, "_duplicate_histogram", ".png", sep=""), width = 10, height = 10)
}

#summary(all_no_na_dups_final[which(all_no_na_dups_final$protein == p & !(all_no_na_dups_final$partner %in% top$partner)),]$delta)


