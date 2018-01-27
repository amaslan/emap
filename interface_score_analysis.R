# plot the various metrics in all df

interface_def = "/Users/annie/emap/alanine_scanning.txt"

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

# for given cluster - protein, boxplot comparison of correlation of correlations whether at interface or not
# determine distribution of correlation of correlations
all_rds = "/Users/annie/emap/all.rds"
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

