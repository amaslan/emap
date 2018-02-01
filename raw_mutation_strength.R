# identify weak mutations
# Generate raw e-map score vs. mutant plots (one plot per complex).
# Only include complexes that have at least one protein with an E-MAP
# score outside of [-3, 2] and WT in [-3, 2].

library(tidyverse)

raw_emap = "/Users/annie/emap/June2016_Gsp1_E-MAP_data.RData"
output_path = "/Users/annie/emap/20180130/"

# load raw e-map data
load(raw_emap)

all <- e.map

# filter: only interested in complexes for which there is at least one protein with 
# min < -3 or max > 2 and WT in [-3, 2]

# remove proteins for which we don't have e-map score
# proteins: 1627 --> 565
# complexes: 408 --> 259
na_removed <- filter(all, !is.na(score))

res <- boxplot(score ~ mutant, data = na_removed)
res_df <- as.data.frame(t(res$stats))
rownames(res_df) <- res$names
colnames(res_df) <- c("lower_whisker", "lower_quartile", "median", "upper_quartile", "upper_whisker")

ggplot(data = res_df) +
  geom_point(mapping = aes(x=lower_whisker, y=upper_whisker)) +
  geom_abline(slope=-1, intercept=0) +
  geom_vline(xintercept=-2, color="red") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(title = "upper and lower whisker for e-map score by mutant; y=-x shown")

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

#strong <- subset(res_df, lower_whisker < -3 | upper_whisker > 2)
strong <- subset(res_df, lower_whisker < -2)

# perform PCA on subset of the data that contains strong muts only
# and exclude WT
strong$mutant <- rownames(strong)
strong_noWT <- subset(strong, !grepl("_WT", mutant))
strong_muts <- rownames(strong_noWT)

# show what -2 cutoff for lower quartile looks like in distribution

ggplot(data=res_df, aes(x=lower_whisker)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(res_df$lower_whisker)), color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=-2), color="red", size=1) +
  labs(x="E-MAP score lower whisker") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(data=res_df, aes(x=upper_whisker)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(res_df$upper_whisker)), color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=1.71), color="red", size=1) +
  labs(x="E-MAP score upper whisker") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

