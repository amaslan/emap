# see if y=-x follows for upper vs. lower whisker for polymerase 
# and histone or if it's unique to our system

library(tidyverse)

# this is the complete ubermap 
pol_emap = "/Users/annie/emap/revised_20180215/polymerase_emap.csv"
hist_emap = "/Users/annie/emap/revised_20180215/histones_pE-MAPs.txt"

####################################################################
# update which 'res <-' line to run based on whether lookint at 
# polymerase or histone data

# polymerase
pol <- read.csv(pol_emap)
pol_no_na <- filter(pol, !is.na(pE.MAP.score))
res <- boxplot(pE.MAP.score ~ RNAPII.Mutant, data = pol_no_na)

# histone
hist <- read.table(hist_emap, 
                       header = TRUE,
                       sep="\t",
                       fill=TRUE,
                       quote="")
hist_no_na <- filter(hist, !is.na(score))
res <- boxplot(score ~ Gene, data = hist_no_na)

####################################################################
# don't need to change anything below here.
# generates boxplots and scatter plots for whiskers & quartiles

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
  labs(title = "upper and lower whisker for e-map score; y=-x shown")

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
