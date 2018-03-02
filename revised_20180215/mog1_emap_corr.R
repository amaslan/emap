# generate simple histograms to motivate why we need correlation 
# of correlations in addition to just raw E-MAP score using
# MOG1 as an example

mog1 <- all[which(all$library_gene_name == 'MOG1'),]
ggplot(data=mog1, aes(x=score)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(mog1$score)), color="blue", linetype="dashed", size=1) +
  labs(x="raw MOG1 E-MAP score") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))

mog1_corr <- sep_no_na[which(sep_no_na$protein == 'MOG1'),]
ggplot(data=mog1_corr, aes(x=value)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(mog1_corr$value)), color="blue", linetype="dashed", size=1) +
  labs(x="MOG1 correlation of correlations") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-1, 1)


### compare to overall

ggplot(data=all, aes(x=score)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(all$score)), color="blue", linetype="dashed", size=1) +
  labs(x="raw all E-MAP score") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-10, 10)

ggplot(data=sep_no_na, aes(x=value)) +
  geom_histogram(aes(y=..density..), color="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") +
  geom_vline(aes(xintercept=mean(sep_no_na$value)), color="blue", linetype="dashed", size=1) +
  labs(x="all correlation of correlations") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  xlim(-1, 1)
