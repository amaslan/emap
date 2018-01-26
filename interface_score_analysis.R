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
interface_short_filtered <- filter(interface_short, score > 0.5 | score < -0.5)

View(interface_short_filtered)
unique(interface_short_filtered$index)
unique(interface_short_filtered$protein)