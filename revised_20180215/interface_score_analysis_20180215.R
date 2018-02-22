library(tidyverse)
library(gplots)
library(RColorBrewer) 
library(rafalib)

# plot the various metrics in all df

interface_def = "/Users/annie/emap/alanine_scanning.txt"
output_path = "/Users/annie/emap/revised_20180215/20180215/"
all_rds = "/Users/annie/emap/revised_20180215/all.rds"

# interface: (1) index	(2) yeast_index	(3) yeast_wt_aa	(4) pdb	(5) protein	(6) structure	(7) pdb_res_num	(8) pdb_wt_aa	(9) score
# e.g.: 1    34          34           T  3a6p    MSN5  MSN5(3a6p)          32       THR 0.000
# just care about index, protein, score
interface <- read.table(interface_def, 
                        header = TRUE,
                        sep="\t",
                        col.names=c("index", "yeast_index", "yeast_wt_aa", "pdb", "protein", "structure", 
                                    "pdb_res_num", "pdb_wt_aa", "score"),
                        fill=TRUE,
                        quote="")
keep_interface = c("index", "protein", "score")
interface_short <- interface[names(interface) %in% keep_interface]
interface_short$index_protein <- paste(interface_short$index, "-", interface_short$protein)

# filter
# interface.score +/-0.5 cutoff
ggplot(data=interface_short, aes(x=score)) +
  geom_histogram(aes(y=..density..), color="black", fill="white", bins=100) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(interface_short$score)), color="blue", linetype="dashed", size=1) +
  labs(x="interface score") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

# this filter results in keeping 14.2% of the data
# residue is at interface with a given protein if it's in interface_short_filtered
interface_short_filtered <- filter(interface_short, score > 0.5 | score < -0.5)
mog1_residues <- c(97, 108, 132, 134, 136, 128, 133, 138)
interface_mog1 <- data.frame(matrix(ncol = 4, nrow = length(mog1_residues)))
colnames(interface_mog1) <- colnames(interface_short_filtered)
interface_mog1$index <- mog1_residues
interface_mog1$protein <- "MOG1"
interface_mog1$score <- NA
interface_mog1$index_protein <- paste(interface_mog1$index, "-", interface_mog1$protein)
interface_short_filtered <- rbind(interface_short_filtered, interface_mog1)

# for given cluster - protein, boxplot comparison of correlation of correlations whether at interface or not
# determine distribution of correlation of correlations
all = readRDS(all_rds)

# remove complex and interface scores info so don't have duplicate correlation of correlation entries
all <- all[,!(names(all) == "Complex")]
all <- all[,!(names(all) == "interface.score")]
all <- all[!duplicated(all),]
# remove NA correlation of correlation values
all$corr.value <- as.numeric(all$corr.value)
all_no_na <- filter(all, !is.na(corr.value))

all_no_na$interface <- all_no_na$residue_protein %in% interface_short_filtered$index_protein

proteins <- c("PSE1", "RNA1", "SRM1", "YRB1", "MOG1")

for (p in proteins) {
  all_no_na_p <- all_no_na %>% filter(protein == p)
  for (c in unique(all_no_na_p$cluster)) {
    ggplot(all_no_na_p[which(all_no_na_p$cluster == c),], aes(x=interface, y=corr.value)) +
      geom_boxplot(colour = "#1F3552", fill = "#4271AE",
                   size = 1) +
      geom_jitter() +
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
    ggsave(filename=paste(output_path, p, "_", c, "_interface_boxplot", ".png", sep=""), width = 10, height = 10)
  }
  }


######################################
# add heatmap of 5 genes vs clusters color by no diff, higher, lower 
# by median or statistical significance value between them

plot_heatmap <- function(df, id, v) {
  mat <- spread(df, cluster, v)
  
  rownames(mat) <- mat[,c(1)]
  mat <- mat[,-1]
  mat <- data.matrix(mat)
  
  hmcol <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
  
  par(mar=c(7,4,4,2)+0.1) 
  png(paste(output_path, id, ".png", sep=""), width=1600, height=1600)
  gplots::heatmap.2(mat, 
                    trace="none",
                    col=hmcol,
                    margins=c(16,14),
                    key.title=NA,
                    main=id)
  dev.off()
}


# protein vs. cluster plot difference of medians
meds <- all_no_na[,c('cluster', 'protein', 'interface', 'corr.value')] %>%
  filter(protein %in% proteins) %>%
  group_by(protein, cluster, interface) %>%
  dplyr::summarize(med = median(corr.value))

m_true <- filter(meds, interface == TRUE)
m_false <- filter(meds, interface == FALSE)
m_diff <- merge(m_true, m_false, by=c('protein', 'cluster'))

# difference in median (interface - non-interface)
m_diff$delta <- m_diff$med.x - m_diff$med.y
  
s <- m_diff[,c('cluster', 'protein', 'delta')]
plot_heatmap(s, 'interface_heatmap_diff_med', 'delta')

s <- m_diff[,c('cluster', 'protein', 'med.x')]
plot_heatmap(s, 'interface_heatmap_interfaceT', 'med.x')

s <- m_diff[,c('cluster', 'protein', 'med.y')]
plot_heatmap(s, 'interface_heatmap_interfaceF', 'med.y')


# statistical significance value instead of just difference in median


# implementation of Mood's test for the median
# dichotomize the data as the pooled median and then Fisher's exact test to see 
# if the binary variable has the same mean in the two samples
median.test <- function(x, y){
  z <- c(x, y)
  g <- rep(1:2, c(length(x), length(y)))
  m <- median(z)
  fisher.test(z < m, g)$p.value
}

df <- data.frame(matrix(ncol=3, nrow=length(proteins)*length(unique(all_no_na_p$cluster))))
colnames(df) <- c('protein', 'cluster', 'pval')
i=1
for (p in proteins) {
  all_no_na_p <- filter(all_no_na, protein == p)
  for (c in unique(all_no_na_p$cluster)) {
    all_no_na_p_c <- all_no_na_p[which(all_no_na_p$cluster == c),]
    interface_true <- all_no_na_p_c[which(all_no_na_p_c$interface == TRUE),]
    interface_false <- all_no_na_p_c[which(all_no_na_p_c$interface == FALSE),]
    df[i,]$pval <- median.test(interface_true$corr.value, interface_false$corr.value)
    df[i,]$protein <- p
    df[i,]$cluster <- c
    i = i+1
  }
}

plot_heatmap(df, 'interface_heatmap_statistical_signficance', 'pval')

