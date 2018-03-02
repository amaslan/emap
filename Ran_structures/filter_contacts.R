# filter contacts to just get interfaces for Ran with binding partners
# interface_residues has final list of residues that are at an interface
# which is used for interface_score_analysis_201802223.R.
# offset of 2 determeine by ran_alignment.R

library(tidyverse)

contacts_txt = "/Users/annie/emap/Ran_structures/contacts.txt"
pdb_filtering_csv = "/Users/annie/emap/Ran_structures/pdb_filtering.csv"

contacts <- read.table(contacts_txt, 
                        header = TRUE,
                        sep="\t",
                        fill=TRUE,
                        quote="")
contacts$complex <- substring(contacts$file_path, 6, 9)
contacts$complex_chain1 <- paste(contacts$complex, "-", contacts$chain1)
contacts$complex_chain2 <- paste(contacts$complex, "-", contacts$chain2)

id <- read.csv(pdb_filtering_csv)
id$complex_chain <- paste(id$complex, "-", id$chain)

# only include the entries where interested in both chains
f <- filter(contacts, contacts$complex_chain1 %in% id$complex_chain & contacts$complex_chain2 %in% id$complex_chain)
f1 <- f[,c('chain1', 'resnum1', 'complex', 'complex_chain1')]
f2 <- f[,c('chain2', 'resnum2', 'complex', 'complex_chain2')]

# merge contacts info with protein id's
f1$complex_chain <- f1$complex_chain1
f2$complex_chain <- f2$complex_chain2
f1_tib <- left_join(as_tibble(f1), as_tibble(id), by='complex_chain')
f2_tib <- left_join(as_tibble(f2), as_tibble(id), by='complex_chain')

# just interested in ran interfaces
f1_ran <- filter(f1_tib, protein == 'RAN')
f2_ran <- filter(f2_tib, protein == 'RAN')

# put all residues in terms of yeast residue numbering
f1_ran_yeast <- filter(f1_ran, organism == 'yeast')
f1_ran_human_canine <- filter(f1_ran, organism == 'human' | organism == 'canine')
f1_ran_human_canine$resnum1 <- f1_ran_human_canine$resnum1 + 2
f2_ran_yeast <- filter(f2_ran, organism == 'yeast')
f2_ran_human_canine <- filter(f2_ran, organism == 'human' | organism == 'canine')
f2_ran_human_canine$resnum2 <- f2_ran_human_canine$resnum2 + 2

# get final list of residue interfaces for ran
interface_residues <- unique(as.vector(rbind(f1_ran_yeast$resnum1, 
                                             f1_ran_human_canine$resnum1,
                                             f2_ran_yeast$resnum2, 
                                             f2_ran_human_canine$resnum2)))

