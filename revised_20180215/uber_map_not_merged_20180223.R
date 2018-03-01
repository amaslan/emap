# see if y=-x follows for upper vs. lower whisker for uber map data
# or if it's unique to our system because it's a switch

library(tidyverse)

# this is the complete ubermap 
uber_map = "/Users/annie/emap/revised_20180215/ubermap_not_merged.RData"

load(uber_map)

# filter to exclude Gene_uniq for which there are fewer than 500 non-NA scores 
# and then UW vs. LW for those

# keep 62% of the data when remove those with score NA
uber_no_na <- filter(not.merged.ubermap.gene.names, !is.na(score))

keep <- uber_no_na %>%
  group_by(query_uniq) %>%
  summarize(total.count = n())

keep_500 <- filter(keep, total.count > 500)
# 4456 genes initially; filtering at >500 for non-NA score count removes 311
length(unique(uber_no_na$query_uniq))
length(unique(keep_500$query_uniq))

ggplot(data=keep, aes(x=total.count)) +
  geom_histogram(aes(y=..density..), color="black", fill="white", bins=100) +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=500), color="red", size=1) +
  labs(x="count of non-NA scores") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

keep_500_genes <- unique(keep_500$query_uniq)

final <- filter(uber_no_na, query_uniq %in% keep_500_genes)

#s <- sample(keep_500_genes, 200)
#sub <- filter(final, query_uniq %in% s)

# just look at point mutations
subPT <- filter(final, !grepl('_damp|_ts', query_uniq, ignore.case = TRUE))
keep_pt <- unique(subPT$query_uniq)
s <- sample(keep_pt, 200)
subPT200 <- filter(subPT, query_uniq %in% s)

res <- boxplot(score ~ query_uniq, data = subPT200)
res_df <- as.data.frame(t(res$stats))
rownames(res_df) <- res$names
colnames(res_df) <- c("lower_whisker", "lower_quartile", "median", "upper_quartile", "upper_whisker")

ggplot(data = res_df, aes(x=lower_whisker, y=upper_whisker)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_abline(slope=-1, intercept=0) +
  ylim(0, 8) +
  xlim(-8, 0) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "upper and lower whisker for e-map score from ubermap data; y=-x shown")


### accessory plots
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

###
# ensure whiskers correlation isn't different for temperature-sensitive
# mutants or knock downs

# check _TS (case insensitive)
subTS <- filter(final, grepl('_ts', query_uniq, ignore.case = TRUE))
length(unique(subTS$query_uniq))

res <- boxplot(score ~ query_uniq, data = subTS)
res_df <- as.data.frame(t(res$stats))
rownames(res_df) <- res$names
colnames(res_df) <- c("lower_whisker", "lower_quartile", "median", "upper_quartile", "upper_whisker")

ggplot(data = res_df, aes(x=lower_whisker, y=upper_whisker)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_abline(slope=-1, intercept=0) +
  ylim(0, 8) +
  xlim(-8, 0) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "TS upper and lower whisker for e-map score from ubermap data; y=-x shown")


# check _DAMP (case insensitive)
subDAMP <- filter(final, grepl('_damp', query_uniq, ignore.case = TRUE))
length(unique(subDAMP$query_uniq))

res <- boxplot(score ~ query_uniq, data = subDAMP)
res_df <- as.data.frame(t(res$stats))
rownames(res_df) <- res$names
colnames(res_df) <- c("lower_whisker", "lower_quartile", "median", "upper_quartile", "upper_whisker")

ggplot(data = res_df, aes(x=lower_whisker, y=upper_whisker)) +
  geom_point() +
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_abline(slope=-1, intercept=0) +
  ylim(0, 8) +
  xlim(-8, 0) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "DAMP upper and lower whisker for e-map score from ubermap data; y=-x shown")
