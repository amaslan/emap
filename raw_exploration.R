# Generate raw e-map score vs. mutant plots (one plot per complex).
# Only include complexes that have at least one protein with an E-MAP
# score outside of [-3, 2] and WT in [-3, 2].

library(tidyverse)

raw_emap = "/Users/annie/emap/June2016_Gsp1_E-MAP_data.RData"
complex_def = "/Users/annie/emap/CYC2008_complex.tab.txt"
output_path = "/Users/annie/emap/20180108/"

# load raw e-map data
load(raw_emap)

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

# merge based on unique identifier ORF
all <- merge(e.map, complex_short, by.x=c("library_ORF"), by.y=c("ORF"), all.y=TRUE)

# filter: only interested in complexes for which there is at least one protein with 
# min < -3 or max > 2 and WT in [-3, 2]

# remove proteins for which we don't have e-map score
# proteins: 1627 --> 565
# complexes: 408 --> 259
na_removed <- filter(all, !is.na(score))

# calculate min e-map score for each protein 
min_score <- na_removed %>% 
  group_by(Name) %>% 
  slice(which.min(score))

# calculate max e-map score for each protein 
max_score <- na_removed %>% 
  group_by(Name) %>% 
  slice(which.max(score))

min_max <- full_join(min_score, max_score, by="Name")

# filter to keep proteins w/ e-map score outside of [-3, 2]
# proteins: 565 --> 452 
# complexes: 259 --> 214 
initial_keep <- min_max[which(min_max$score.x < -3 | min_max$score.y > 2),]
protein_keep = initial_keep$Name

# now need exclude proteins if WT outside of [-3, 2]
# proteins: 452 --> 442
# complexes: 214 --> 220 (some proteins in multiple complexes)
wt <- filter(all, grepl("-NAT", mutant))
wt_keep <- filter(wt, Name %in% protein_keep)
final <- wt_keep[which(wt_keep$score >= -3 & wt_keep$score <= 2),]

# save all the plots as pngs with one per complex (220 after all filtering)
i=0
for (c in unique(final$Complex)) {
  ggplot(data = all[which(all$Complex == c),]) +
          geom_point(mapping = aes(x=reorder(mutant, score), y=score, color=Name)) +
          labs(x="mutant", y="E-MAP score", title=c) +
          geom_hline(yintercept = 0) +
          theme(axis.text.x = element_text(angle = 90, hjust=1), panel.background = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(filename=paste(output_path, i, "_raw_emap", ".png", sep=""), width = 10, height = 10)
  i=i+1
}

