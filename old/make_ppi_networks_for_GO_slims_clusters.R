### this code is based on https://www.r-bloggers.com/network-visualization-part-5-cytoscape-an-update-rcy3/

#### this script should combine the clustering data (GO-slims pearson complete clusters I use for corr_of_corr E-MAP analysis)
#### with the biogrid data on protein-protein interactions in yeast and makes cytoscape networks

#### get all proteins from the cluster and draw all the ppi edges between them
#### then add a "background network" that will contain all the immediate partners of all the cluster proteins (add all the edges too)
library(igraph)
library(graph)
library(RCy3)
library(RJSONIO)
options(stringsAsFactors = F)

ppi_df <- read.delim("/Users/annie/emap/Scerevisiae_biogrid_preprocessed_simple_high_confidence_heteromers.txt", head = T)
#clusters_df <- read.delim("2016-11-02_GO_slims_pearson_complete_clusters.txt", head = T)
clusters_df <- sep

clusters <- as.character(unique(clusters_df$cluster))
for ( i in seq_along(clusters) ) {
  by_cluster_df <- clusters_df[clusters_df$cluster == clusters[i],]
  by_cluster_ORFs <- as.character(unique(by_cluster_df$ORF))
  by_cluster_genes <- as.character(unique(by_cluster_df$gene_name))
  #### get all high confidence ppi involving genes from this cluster
  by_cluster_ppi <- ppi_df[ppi_df$ORF1 %in% by_cluster_ORFs | ppi_df$ORF2 %in% by_cluster_ORFs,]
  ### only the high confidence ppi between cluster members
  #by_cluster_ppi_core <- by_cluster_ppi[
   # by_cluster_ppi$ORF1 %in% by_cluster_ORFs & by_cluster_ppi$ORF2 %in% by_cluster_ORFs,]
  by_cluster_ppi <- by_cluster_ppi[, c(3,4,6)]
  # make it into a graph object
  by_cluster_ppi.network <- graph.data.frame(by_cluster_ppi, directed = F)
  ##### make node attributes - is it in the cluster or not
  by_cluster_ppi.network <- set_vertex_attr( by_cluster_ppi.network, 
      'label', index = V(by_cluster_ppi.network)$name %in% by_cluster_genes, value = "cluster")
  by_cluster_ppi.network <- set_vertex_attr( by_cluster_ppi.network, 
      'label', index = (! V(by_cluster_ppi.network)$name %in% by_cluster_genes), value = "outside")
  for (j in 1:nrow(by_cluster_ppi)) {
    igraph::E(by_cluster_ppi.network)[
        as.character(by_cluster_ppi$gene1) %--% as.character(by_cluster_ppi$gene2)]$weight <- 
      as.numeric(by_cluster_ppi$unique_experiment_count)
  }  ### %--% is a special operator from igraph
  # make it into a graphNEL object
  by_cluster_ppi.cyt <- igraph:::as_graphnel(by_cluster_ppi.network)
  # the same with core
  #by_cluster_ppi_core <- by_cluster_ppi_core[, c(3,4,6)]
  #by_cluster_ppi_core.network <- graph.data.frame(by_cluster_ppi_core, directed = F)
  #by_cluster_ppi_core.cyt <- igraph::as_graphnel(by_cluster_ppi_core.network)
  by_cluster_ppi.cyt <- RCy3::initNodeAttribute(by_cluster_ppi.cyt, 'label', 'char', 'none') 
  #by_cluster_ppi.cyt <- RCy3::initNodeAttribute(by_cluster_ppi.cyt, 'name', 'char', 'GENE')
  by_cluster_ppi.cyt <- RCy3::initEdgeAttribute (by_cluster_ppi.cyt, "weight", 'integer', 0)
  by_cluster_ppi.cyt <- RCy3::initEdgeAttribute (by_cluster_ppi.cyt, "unique_experiment_count", 'integer', 2)
  gDCW <- RCy3::CytoscapeWindow(clusters[i], graph = by_cluster_ppi.cyt, overwriteWindow = TRUE)
  RCy3::displayGraph(gDCW)
  RCy3::setNodeAttributesDirect(gDCW, 'label', 'char', igraph::V(by_cluster_ppi.network)$name, igraph::V(by_cluster_ppi.network)$label)
  RCy3::setNodeAttributesDirect(gDCW, 'name', 'char', igraph::V(by_cluster_ppi.network)$name, igraph::V(by_cluster_ppi.network)$name)
  RCy3::setEdgeAttributesDirect(gDCW, 'weight', 'integer', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(by_cluster_ppi.cyt, as.character(by_cluster_ppi$gene1), as.character(by_cluster_ppi$gene2), 'weight'))
  RCy3::setEdgeAttributesDirect(gDCW, 'unique_experiment_count', 'numeric', as.character (RCy3::cy2.edge.names (gDCW@graph)), graph::edgeData(by_cluster_ppi.cyt, as.character(by_cluster_ppi$gene1), as.character(by_cluster_ppi$gene2), 'unique_experiment_count'))
  cy <- RCy3::CytoscapeConnection()
  hlp <-RCy3::getLayoutNames(cy)
  getLayoutPropertyNames(cy, hlp[7])
  RCy3::setLayoutProperties (gDCW, hlp[7], list (numIterations = 1000))
  RCy3::layoutNetwork(gDCW, hlp[7])
  RCy3::setNodeColorRule(gDCW, 'label', c("cluster", "outside"), c('#FF7F50', '#F5DEB3'), mode = 'lookup')
  RCy3::setNodeSizeRule(gDCW, 'label', c("cluster", "outside"), c(60, 40), mode = 'lookup')
  # And edges:
  RCy3::setEdgeLineWidthRule(gDCW, 'weight', by_cluster_ppi$unique_experiment_count, by_cluster_ppi$unique_experiment_count)
  RCy3::setDefaultBackgroundColor(gDCW, '#D3D3D3')
  RCy3::setDefaultNodeColor(gDCW, '#87CEFA')
  RCy3::setDefaultNodeLabelColor(gDCW, '#000000')
  file.name <- paste(clusters[i], "_network")
  saveImage(gDCW, file.name, image.type = 'pdf')
}





