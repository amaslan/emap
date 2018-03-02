# whisker correlation plots for subset of ubermap data that overlaps
# the library our mutants were screened against
# NB. script takes subset of 200 genes at a time so R session doesn't 
# time out

library(tidyverse)

# this is the part of the ubermap that overlaps with the library our mutants were screened against
uber_map = "/Users/annie/emap/revised_20180215/preprocessed_ubermap_ubergenes_only.txt"

uber <- read.table(uber_map, 
                        header = TRUE,
                        sep="\t",
                        fill=TRUE,
                        quote="")

# filter to exclude Gene_uniq for which there are fewer than 500 non-NA scores 
# and then UW vs. LW for those

# keep 62% of the data when remove those with score NA
uber_no_na <- filter(uber, !is.na(score))

keep <- uber_no_na %>%
  group_by(Gene_uniq) %>%
  summarize(total.count = n())

keep_500 <- filter(keep, total.count > 500)
# 4456 genes initially; filtering at >500 for non-NA score count removes 311
length(unique(uber_no_na$Gene_uniq))
length(unique(keep_500$Gene_uniq))

ggplot(data=keep, aes(x=total.count)) +
  geom_histogram(aes(y=..density..), color="black", fill="white", bins=100) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=500), color="red", size=1) +
  labs(x="count of non-NA scores") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

keep_500_genes <- unique(keep_500$Gene_uniq)

# this filtering results in keeping 97% of the data
final <- filter(uber_no_na, Gene_uniq %in% keep_500_genes)

# taking 200 randomly at a time so R session doesn't timeout
s <- sample(keep_500_genes, 200)
sub <- filter(final, Gene_uniq %in% s)
  
res <- boxplot(score ~ Gene_uniq, data = sub)
res_df <- as.data.frame(t(res$stats))
rownames(res_df) <- res$names
colnames(res_df) <- c("lower_whisker", "lower_quartile", "median", "upper_quartile", "upper_whisker")

ggplot(data = res_df, aes(x=lower_whisker, y=upper_whisker)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_abline(slope=-1, intercept=0) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "upper and lower whisker for e-map score from ubermap data; y=-x shown")

ggplot(data = res_df) +
  geom_point(mapping = aes(x=lower_quartile, y=upper_quartile)) +
  geom_abline(slope=-1, intercept=0) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "upper and lower quartile for e-map score by mutant; y=-x shown")

ggplot(data = res_df) +
  geom_point(mapping = aes(x=upper_whisker, y=upper_quartile)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "upper quartile and upper whisker for e-map score by mutant")

ggplot(data = res_df) +
  geom_point(mapping = aes(x=lower_whisker, y=lower_quartile)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "lower quartile and lower whisker for e-map score by mutant")



