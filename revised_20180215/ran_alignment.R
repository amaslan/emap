# align yeast GSP1 / human RAN sequences for translate residues
# http://www.uniprot.org/uniprot/P62826
# starting at residue 8 for human and 10 for yeast, offset is just +2

library(Biostrings)

yeast = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL'
human = 'MAAQGEPQVQFKLVLVGDGGTGKTTFVKRHLTGEFEKKYVATLGVEVHPLVFHTNRGPIKFNVWDTAGQEKFGGLRDGYYIQAQCAIIMFDVTSRVTYKNVPNWHRDLVRVCENIPIVLCGNKVDIKDRKVKAKSIVFHRKKNLQYYDISAKSNYNFEKPFLWLARKLIGDPNLEFVAMPALAPPEVVMDPALAAQYEHDLEVAQTTALPDEDDDL'

p <- pairwiseAlignment(pattern = c(yeast), subject = human)

frac_match = nmatch(p)/nchar(human)
perc_sequence_identity = pid(p)

# yeast residue numbers are 2 larger than human starting at human 8 / yeast 10
p@pattern@mismatch@unlistData
p@subject@mismatch@unlistData

p@pattern@mismatch@unlistData - p@subject@mismatch@unlistData

