# create contacts.txt
# remove chains that aren't of interest
# change chain for GTPs
# add gene names

library(bio3d)
library(tidyverse)
setwd("/Users/annie/emap/Ran_structures")
pdb_dir <- "pdbs"
files_list <- list.files(path = pdb_dir)
file_paths <- file.path(pdb_dir, files_list)
contacts_table <- data.frame()
for (i in seq_along(file_paths)) {
  file_path <- file_paths[i]
  pdb <- read.pdb(file_path, rm.alt = FALSE)
  coord.temp <- as_tibble(pdb$atom[, c("chain", "resno", "resid", "elety", "x", "y", "z")]) %>%
    filter(resid != "HOH")
  chains <- unique(coord.temp$chain)
  coords <- list()
  for (ch in seq_along(chains)) {
    chain_id <- chains[ch]
    coords[[chain_id]] <- coord.temp %>% 
      filter(chain == chain_id) %>% 
      mutate("unique" = str_c(chain, resno, elety, sep = "_"))
  }
  pairs_of_chains <- combn(chains, 2)
  for (p in seq_along(pairs_of_chains[1,])) {
    chain1 <- pairs_of_chains[1, p]
    chain2 <- pairs_of_chains[2, p]
    matrix1 <- as.matrix(coords[[chain1]][, c("x", "y", "z")])
    rownames(matrix1) <- coords[[chain1]]$unique
    matrix2 <- as.matrix(coords[[chain2]][, c("x", "y", "z")])
    rownames(matrix2) <- coords[[chain2]]$unique
    distances <- dist.xyz(matrix1, matrix2)
    colnames(distances) <- rownames(matrix2)
    distances_df <- data.frame(distances)
    distances_df <- cbind(distances_df, "unique_chain1" = unique(as.character(rownames(matrix1))))
    contacts_df <- as_tibble(distances_df) %>% 
      gather(unique_chain2, dist, -unique_chain1) %>%
      filter(dist < 4) %>%
      separate(col = unique_chain1, into = c("chain1", "resnum1", "atom1"), convert = T) %>%
      separate(col = unique_chain2, into = c("chain2", "resnum2", "atom2"), convert = T) %>%
      group_by(resnum1, resnum2) %>%
      dplyr::summarise("number_of_contacts" = n()) %>%
      arrange(resnum1, resnum2)
    if (nrow(contacts_df) > 0) {
      contacts_table <- rbind(contacts_table, as_tibble(data.frame(file_path, chain1, chain2, contacts_df)))
    }
  }
}
contacts_table <- as_tibble(contacts_table)
write_delim(contacts_table, "contacts.txt", delim = "\t")

# filter to look at residues that are 