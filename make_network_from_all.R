# data and settings for cytoscape graph

library(tidyverse)

# determine distribution of correlation of correlations
all_rds = "/Users/annie/emap/all.rds"
output_path = "/Users/annie/emap/20180116/"

all = readRDS(all_rds)

# remove complex info so don't have duplicate correlation of correlation entries
all <- all[,!(names(all) == "Complex")]
all <- all[!duplicated(all),]

all$corr.value <- as.numeric(all$corr.value)
all_no_na <- filter(all, !is.na(corr.value))

# plot histogram of correlation of correlation values
ggplot(data=all_no_na, aes(x=corr.value)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(all_no_na$corr.value)), color="blue", linetype="dashed", size=1) +
  labs(x="correlation of correlations") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

# cutoff for correlation of correlation values  
  
# p = pnorm(result_m$value, mean = mean(result_m$value), sd = sd(result_m$value)/sqrt(length(result_m$value)))
# pBH = p.adjust(p, method="BH")
# pBonferroni = p.adjust(p, method="bonferroni")
# 
# plot(p[1:100], pBonferroni[1:100])
# lines(lowess(p[1:100],pBonferroni[1:100]), col="blue")
# points(p[1:100], pBH[1:100])
# points(p[1:100], p[1:100])

# plot histogram of e.map values
ggplot(data=all_no_na, aes(x=raw.emap)) +
  geom_histogram(aes(y=..density..), color="black", fill="white", bins = 100) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(all_no_na$raw.emap, na.rm=TRUE)), color="blue", linetype="dashed", size=1) +
  labs(x="raw E-MAP score") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

# plot histogram of interface scores
ggplot(data=all_no_na, aes(x=interface.score)) +
  geom_histogram(aes(y=..density..), color="black", fill="white", bins = 100) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(all_no_na$interface.score, na.rm=TRUE)), color="blue", linetype="dashed", size=1) +
  labs(x="abs val interface score") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

######################################
### this code is based on https://www.r-bloggers.com/network-visualization-part-5-cytoscape-an-update-rcy3/ & Tina Perica's code

### this script uses clustering data (GO-slims pearson complete clusters I use for corr_of_corr E-MAP analysis) to create
### cytoscape networks of protein-mutant interactions in yeast

### get all proteins from the cluster and draw all edges to proteins/mutants weighted by corr_of_corr

library(igraph)
library(graph)
library(RCy3)
library(RJSONIO)

# summarize duplicates for mutant - protein
# dup_genes contains the 18 genes for which there are multiple experiments
dups <- aggregate(all_no_na$gene_mut, list(all_no_na$cluster, all_no_na$protein), length)
dups <- dups[which(dups$x > length(unique(all_no_na$mutant))),]
dup_genes <- unique(dups$Group.2)


# whether considering experiments with same gene-mutant combo separately or together
clusters_df <- all_no_na[which(!(all_no_na$protein %in% dup_genes)),]

# will use statistical test to determine which edges to include; for testing abs val corr of corr > 0.6
clusters_df <- clusters_df[which(abs(clusters_df$corr.value) > 0.6),]

clusters <- as.character(unique(clusters_df$cluster))

i=1
for ( i in seq_along(clusters) ) {
  by_cluster_df <- clusters_df[clusters_df$cluster == clusters[i],]
  by_cluster_ORFs <- as.character(unique(by_cluster_df$ORF))
  by_cluster_genes <- as.character(unique(by_cluster_df$protein))
  
  # unlike Tina's script, here we consider all genes from this cluster (not just high confidence ppi)
  # required data format is dataframe with 3 variables; variables 1&2 correspond to interactions; variable 3 is weight of interaction
  by_cluster_F <- by_cluster_df[,c('mutant', 'protein', 'corr.value')]
  
  # make it into a graph object
  by_cluster_F.network <- graph.data.frame(by_cluster_F, directed = F)
  
  ##### make node attributes - is it mutant or gene
  by_cluster_F.network <- set_vertex_attr( by_cluster_F.network, 
                                             'label', index = V(by_cluster_F.network)$name %in% by_cluster_genes, value = "gene")
  by_cluster_F.network <- set_vertex_attr( by_cluster_F.network, 
                                             'label', index = (! V(by_cluster_F.network)$name %in% by_cluster_genes), value = "mutant")
  
  # asign edge weights based on correlation of correlations
  for (j in 1:nrow(by_cluster_F)) {
    igraph::E(by_cluster_F.network)[
      as.character(by_cluster_F$mutant) %--% as.character(by_cluster_F$protein)]$weight <- 
      as.numeric(by_cluster_F$corr.value)
  }  ### %--% is a special operator from igraph
  
  # make it into a graphNEL object
  by_cluster_F.cyt <- igraph:::as_graphnel(by_cluster_F.network)
  
  by_cluster_F.cyt <- RCy3::initNodeAttribute(by_cluster_F.cyt, 'label', 'char', 'none') 
  by_cluster_F.cyt <- RCy3::initNodeAttribute(by_cluster_F.cyt, 'name', 'char', 'GENE/MUT')
  by_cluster_F.cyt <- RCy3::initEdgeAttribute (by_cluster_F.cyt, "weight", 'integer', 0)
  by_cluster_F.cyt <- RCy3::initEdgeAttribute (by_cluster_F.cyt, "corr.value", 'integer', 2)
  
  gDCW <- RCy3::CytoscapeWindow(clusters[i], graph = by_cluster_F.cyt, overwriteWindow = TRUE)
  
  RCy3::displayGraph(gDCW)
  
  RCy3::setNodeAttributesDirect(gDCW, 'label', 'char', igraph::V(by_cluster_F.network)$name, igraph::V(by_cluster_F.network)$label)
  RCy3::setNodeAttributesDirect(gDCW, 'name', 'char', igraph::V(by_cluster_F.network)$name, igraph::V(by_cluster_F.network)$name)
  RCy3::setEdgeAttributesDirect(gDCW, 'weight', 'integer', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(by_cluster_F.cyt, as.character(by_cluster_F$mutant), as.character(by_cluster_F$protein), 'weight'))
  RCy3::setEdgeAttributesDirect(gDCW, 'corr.value', 'numeric', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(by_cluster_F.cyt, as.character(by_cluster_F$mutant), as.character(by_cluster_F$protein), 'corr.value'))
  
  cy <- RCy3::CytoscapeConnection()
  hlp <-RCy3::getLayoutNames(cy)
  getLayoutPropertyNames(cy, hlp[7])
  RCy3::setLayoutProperties (gDCW, hlp[7], list (numIterations = 1000))
  RCy3::layoutNetwork(gDCW, hlp[7])
  RCy3::setNodeColorRule(gDCW, 'label', c("protein", "mutant"), c('#FF7F50', '#F5DEB3'), mode = 'lookup')
  RCy3::setNodeSizeRule(gDCW, 'label', c("protein", "mutant"), c(60, 40), mode = 'lookup')
  # And edges:
  #RCy3::setEdgeLineWidthRule(gDCW, 'weight', by_cluster_ppi$unique_experiment_count, by_cluster_ppi$unique_experiment_count)
  RCy3::setDefaultBackgroundColor(gDCW, '#D3D3D3')
  RCy3::setDefaultNodeColor(gDCW, '#87CEFA')
  RCy3::setDefaultNodeLabelColor(gDCW, '#000000')
  file.name <- paste(clusters[i], "_network")
  saveImage(gDCW, file.name, image.type = 'pdf')
}













