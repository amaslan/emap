# plot raw E-MAP score vs. MDS distance or number of edges between genes
# do we produce a similar plot to Collins, et al.?

library(tidyverse)

raw_emap = "/Users/annie/emap/June2016_Gsp1_E-MAP_data.RData"
corr_of_corr = "/Users/annie/emap/corr_of_corr.RData"

complex_def = "/Users/annie/emap/CYC2008_complex.tab.txt"
interface_def = "/Users/annie/emap/alanine_scanning.txt"

method = "corr_of_corr_no_na_no_mut"

########################################
### create giant dataframe with all info: 
### raw e-map, corr of corr, complex, interface

### read in all data
# e.map: (1) mutant (2) library_ORF (3) library_gene_name (4) score
# e.g.: 1         T137G     YAL002W              VPS8  0.117158
load(raw_emap)

# result$corr_of_corr: (1) method (2) value (3) mutant (4) partner (5) cluster_number (6) cluster
# e.g.: 1      corr_of_corr_all  -0.423138711324343  K129F CIN1 - YOR349W             10 Golgi and ER_GO_2
load(corr_of_corr)
result_m <- result$corr_of_corr[which(result$corr_of_corr$method == method),]
sep <- result_m %>% 
  separate(partner, into=c("protein", "ORF"), sep=" - ", remove=FALSE)

# complex: (1) ORF	(2) Name	(3) Complex	(4) PubMed_id	(5) Method	(6) GO_id	(7) GO_term	(8) Jaccard_Index
# e.g.: 1 YKR068C   BET3 TRAPP complex  10727015 "Affinity Capture-Western,Affinity Capture-MS" GO:0030008 TRAPP complex             1
# read in file that defines complexes
complex <- read.table(complex_def, 
                      header = TRUE,
                      sep="\t",
                      col.names=c("ORF", "Name", "Complex", "PubMed_id", "Method", "GO_id", 
                                  "GO_term", "Jaccard_Index"),
                      fill=TRUE,
                      quote="")

# just care about ORF, Name, Complex from complex table
keep = c("ORF", "Name", "Complex")
complex_short <- complex[names(complex) %in% keep]

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

### merge and save all data
# merge raw e.map and complex_short based on unique identifier ORF
all <- merge(e.map, complex_short, by.x=c("library_ORF"), by.y=c("ORF"), all=TRUE)

# merge all with sep from corr_of_corr based on unique mutant-gene combo
all$gene_mut <- paste(all$Name, "-", all$mutant)
sep$gene_mut <- paste(sep$protein, "-", sep$mutant)
all <- merge(all, sep, by="gene_mut", all=TRUE)

# merge all with interface based on residue number + protein
interface_short$residue_protein <- paste(interface_short$index, "-", interface_short$protein)
by_residue_protein <- group_by(interface_short, residue_protein)
interface_clean <- summarize(by_residue_protein, score_max = max(abs(score)))

all$residue <- parse_number(all$mutant.y)
all$residue_protein <- paste(all$residue, "-", all$protein)
all <- merge(all, interface_clean, by="residue_protein", all=TRUE)

# cleanup and save dataframe
drop <- c("library_ORF", "mutant.x", "library_gene_name", "Name")
all <- all[,!(names(all) %in% drop)]
names(all)[names(all) == "mutant.y"] <- "mutant"
names(all)[names(all) == "score"] <- "raw.emap"
names(all)[names(all) == "score_max"] <- "interface.score"
names(all)[names(all) == "value"] <- "corr.value"

all$gene_mut_cluster <- paste(all$gene_mut, "-", all$cluster)
all <- all[!duplicated(all),]
saveRDS(all, file="/Users/annie/emap/all.rds")

########################################

# # filter
# # interface.score +/-0.5 cutoff
# ggplot(data=all, aes(x=interface.score)) +
#   geom_histogram(aes(y=..density..), color="black", fill="white") +
#   geom_density(alpha=.2, fill="#FF6666") +
#   geom_vline(aes(xintercept=mean(all$interface.score)), color="blue", linetype="dashed", size=1) +
#   labs(x="interface score") +
#   theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# # this filter results in keeping 10.9% of the data
# all_interface <- filter(all, interface.score > 0.5 | interface.score < -0.5)
# 
# # filter for raw e.map scores > 2 or < -3
# all_raw <- filter(all, raw.emap > 2 | raw.emap < -3)
#   
# plot(all_interface$interface.score, all_interface$corr.value)
# plot(all$corr.value, all$raw.emap)
# 
# ggplot(data = all_raw) + 
#   #geom_point(mapping = aes(x = displ, y = hwy)) +
#   geom_smooth(mapping = aes(x = corr.value, y = raw.emap), na.rm = TRUE)
# 
# plot(all_raw$corr.value, all_raw$raw.emap)


