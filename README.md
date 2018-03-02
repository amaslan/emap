https://github.com/amaslan/emap
## summary of most useful scripts
### emap/revised_20180215/

- cluster_size.R - calculate number of genes in each cluster for real and random clustering
- corr_of_corr_exploration_20180215.R - generate heatmaps: protein vs. mutant (1 per cluster), cluster vs. mutant (1 per protein), protein-cluster vs. mutant
- dup_gene_analysis_20180215.R - look at genes with multiple experiments (28) and determine distribution of difference in correlation of correlations, etc. → used to determine that we will consider genes separately
- interface_score_analysis_20180215.R - interface analysis using alanine scanning data → generate boxplots for correlation of correlations values at interfaces vs. not and heatmaps of difference in median and significance of median differences for interface vs. not
- interface_score_analysis_20180223.R - analysis of correlation of correlations values separate by interface vs. not from pdb - must run /emap/Ran_structures/filter_contacts.R first to get interfaces_residues vector for this script
- make_all_df_20180215.R - create single dataframe with: (1) raw e-map, (2) correlation of correlations, (3) complex, (4) interface from alanine scanning
- mog1_emap_corr.R - generate MOG1 raw E-MAP and correlation of correlations histograms (and do for all genes as well) to motivate why need correlation of correlations
- mut_clustering_20180228.R - cluster mutations → generate dendrograms and do MDS for mutant pairs (first pairwise distances for mutant given vector of correlation of correlations with proteins; then distances btw pairs across 34 clusters)
- pca_20180227.R - PCA of correlation of correlations values to see if we can distinguish sub-GO categories from hierarchical clustering & if we can distinguish GO categories - need to specify which subset of mutations and which subset of genes want to consider in first section of the script
- polymerase_histone_20180227.R - produce whisker correlation plots for polymerase and histone data to see if we see y=-x trend
- ran_alignment.R - small script to get offset for Gsp1/Ran residues for yeast / human / canine
- uber_map_20180221.R - produce whisker correlation plots for subset of ubermap data that overlaps the library our mutants were screened against (NB script currently takes random subset of 200 genes at a time so  R session doesn’t time out)
- uber_map_not_merged_20180223.R - produce whisker correlation plots for all ubermap data (NB script currently takes random subset of 200 genes at a time so  R session doesn’t time out); produce separate plots for gene knockout, gene knock down, and temperature sensitive mutants
### emap/
- interface_score_analysis.R - produce boxplots of correlation of correlations values separated by interface vs. not for each cluster. Also makes heatmaps of median difference in correlation of correlations by at interface or not and significance of difference. Interface defined from alanine scanning.
- make_network_from_all_final.R - create cytoscape graph; older analysis so was original correlation of correlations values, excludes duplicate genes, and has arbitrary correlation of correlations value cutoff; calculates betweenness centrality and degree & summary statistics
- make_network_from_all_final_strong.R - same as above but for strong mutants only and 5 genes (RNA1, SRM1, MOG1, YRB1, PSE1) only
- raw_exploration.R - generate titration curves of raw e-map score vs. mutant plots (1 plot per complex). Only include complexes that have at least one protein with an E-MAP score outside of [-3, 2] and WT in [-3, 2]
- raw_mutation_strength.R - generate whisker correlation plot and analysis to determine threshold for strong muts to consider in downstream analysis

### emap/Ran_structures
- filter_contacts.R - filter contacts to just get interfaces for Ran with binding partners