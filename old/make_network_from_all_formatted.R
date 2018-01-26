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
library(plyr)

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
  by_cluster_genes <- as.character(unique(by_cluster_df$protein))
  
  # required data format is dataframe with 3 variables; variables 1&2 correspond to interactions; variable 3 is weight of interaction
  by_cluster_short <- by_cluster_df[,c('mutant', 'protein', 'corr.value')]
  
  # make it into a graph object - use simplify to ensure that there are no duplicated edges or self loops
  gD <- simplify(graph.data.frame(by_cluster_short, directed = FALSE))
  
  ######################################################################
  ### calculate node properties
  # Calculate degree for all nodes
  degAll <- igraph::degree(gD, v = V(gD), mode = "all")
  
  # Calculate betweenness for all nodes
  betAll <- betweenness(gD, v = V(gD), directed = FALSE) / (((vcount(gD) - 1) * (vcount(gD)-2)) / 2)
  betAll.norm <- (betAll - min(betAll))/(max(betAll) - min(betAll))
  rm(betAll)
  
  # Calculate Dice similarities between all pairs of nodes
  dsAll <- similarity.dice(gD, vids = V(gD), mode = "all")
  
  ######################################################################
  ### Add new node/edge attributes based on the calculated node properties/similarities
  
  gD <- set.vertex.attribute(gD, "degree", index = V(gD), value = degAll)
  gD <- set.vertex.attribute(gD, "betweenness", index = V(gD), value = betAll.norm)
  gD <- set.edge.attribute(gD, "weight", index = E(gD), value = 0)
  gD <- set.edge.attribute(gD, "similarity", index = E(gD), value = 0)
  gD <- set.vertex.attribute(gD, "label", index = V(gD)$name %in% by_cluster_genes, value = "gene")
  gD <- set.vertex.attribute(gD, "label", index = (! V(gD)$name %in% by_cluster_genes), value = "mutant")
  
  F1 <- function(x) {data.frame(sim = dsAll[which(V(gD)$name == as.character(x$mutant)), which(V(gD)$name == as.character(x$protein))])}
  dataSet.ext <- ddply(by_cluster_short, .variables=c("mutant", "protein", "corr.value"), function(x) data.frame(F1(x)))
  
  # The order of interactions in gD is not the same as it is in dataSet or as it is in the edge list,
  # and for that reason these values cannot be assigned directly
  ### %--% is a special operator from igraph
  E(gD)[as.character(dataSet.ext$mutant) %--% as.character(dataSet.ext$mutant)]$weight <- as.numeric(dataSet.ext$corr.value)
  E(gD)[as.character(dataSet.ext$mutant) %--% as.character(dataSet.ext$mutant)]$similarity <- as.numeric(dataSet.ext$sim)
  
  ##################################################
  ### Print network in Cytoscape
  # make it into a graphNEL object
  gD.cyt <- igraph.to.graphNEL(gD)
  
  # We have to create attributes for graphNEL
  # We'll keep the same name, so the values are passed from igraph
  
  gD.cyt <- initNodeAttribute(gD.cyt, 'degree', 'numeric', 0) 
  gD.cyt <- initNodeAttribute(gD.cyt, 'betweenness', 'numeric', 0) 
  gD.cyt <- initEdgeAttribute(gD.cyt, "weight", 'integer', 0)
  gD.cyt <- initEdgeAttribute(gD.cyt, "similarity", 'numeric', 0)
  gD.cyt <- initNodeAttribute(gD.cyt, 'label', 'char', 'none') 
  
  # Now we can create a new graph window in cytoscape and display graph with default color/size scheme
  gDCW <- CytoscapeWindow(clusters[i], graph = gD.cyt, overwriteWindow = TRUE)
  displayGraph(gDCW)
  
  setNodeAttributesDirect(gDCW, 'label', 'char', igraph::V(gD)$name, igraph::V(gD)$label)
  #setNodeAttributesDirect(gDCW, 'name', 'char', igraph::V(by_cluster_F.network)$name, igraph::V(by_cluster_F.network)$name)
  setEdgeAttributesDirect(gDCW, 'weight', 'integer', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(gD.cyt, as.character(gD$mutant), as.character(gD$protein), 'weight'))
  setEdgeAttributesDirect(gDCW, 'corr.value', 'numeric', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(gD.cyt, as.character(gD$mutant), as.character(gD$protein), 'corr.value'))
  
  
  # Apply values to some of the properties and plot the layout
  cy <- CytoscapeConnection()
  hlp <-getLayoutNames(cy)
  setLayoutProperties(gDCW, hlp[7], list(nIterations = 1000))
  layoutNetwork(gDCW, hlp[7])
  
  setNodeColorRule(gDCW, 'label', c("protein", "mutant"), c('#FF7F50', '#F5DEB3'), mode = 'lookup')
  setNodeSizeRule(gDCW, 'label', c("protein", "mutant"), c(60, 40), mode = 'lookup')
  
  # Now, we can define our own default color/size scheme
  # setDefaultBackgroundColor(gDCW, '#FFFFFF')
  # setDefaultEdgeColor(gDCW, '#CDC9C9')
  # setDefaultEdgeLineWidth(gDCW, 4)
  # setDefaultNodeBorderColor(gDCW, '#000000')
  # setDefaultNodeBorderWidth(gDCW, 3)
  # setDefaultNodeShape(gDCW, 'ellipse')
  # setDefaultNodeColor(gDCW, '#87CEFA')
  # setDefaultNodeSize(gDCW, 60)
  # setDefaultNodeFontSize(gDCW, 20)
  # setDefaultNodeLabelColor(gDCW, '#000000')
  
  
  # Finally, we can define rules for node colors, node sizes, and edge colors
  setNodeColorRule(gDCW, 'degree', c(min(degAll), mean(degAll), max(degAll)), c('#F5DEB3', '#FFA500', '#FF7F50', '#FF4500', '#FF0000'), mode = 'interpolate')
  setNodeSizeRule(gDCW, 'betweenness', c(min(betAll.norm), mean(betAll.norm), max(betAll.norm)), c(30, 45, 60, 80, 100), mode = 'interpolate')
  #setEdgeColorRule(gDCW, 'weight', c(min(as.numeric(dataSet.ext$corr.value)), mean(as.numeric(dataSet.ext$corr.value)), max(as.numeric(dataSet.ext$corr.value))), c('#FFFF00', '#00FFFF', '#00FF7F', '#228B22', '#006400'), mode='interpolate')
  redraw(gDCW)
  
  setEdgeLineWidthRule(gDCW, 'weight', gD$corr.value, gD$corr.value)
  
  
  setNodeAttributesDirect(gDCW, 'label', 'char', igraph::V(gD)$name, igraph::V(gD)$label)
  setNodeAttributesDirect(gDCW, 'degree', 'numeric', igraph::V(gD)$name, igraph::V(gD)$degree)
  setNodeAttributesDirect(gDCW, 'betweenness', 'numeric', igraph::V(gD)$name, igraph::V(gD)$betweenness)
  
  setEdgeAttributesDirect(gDCW, 'weight', 'integer', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(gD.cyt, as.character(gD$mutant), as.character(gD$protein), 'weight'))
  setEdgeAttributesDirect(gDCW, 'similarity', 'numeric', igraph::V(gD)$name, igraph::V(gD)$similarity)
  
  
  RCy3::setEdgeAttributesDirect(gDCW, 'weight', 'integer', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(gD.cyt, as.character(by_cluster_F$mutant), as.character(by_cluster_F$protein), 'weight'))
  RCy3::setEdgeAttributesDirect(gDCW, 'corr.value', 'numeric', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(gD.cyt, as.character(by_cluster_F$mutant), as.character(by_cluster_F$protein), 'corr.value'))
  
  
  file.name <- paste(clusters[i], "_network")
  saveImage(gDCW, file.name, image.type = 'pdf')
}













