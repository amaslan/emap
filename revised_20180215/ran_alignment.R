# align yeast GSP1 / human / canine RAN sequences for translate residues
# http://www.uniprot.org/uniprot/P62826
# starting at residue 8 for human/canine and 10 for yeast, offset is just +2

library(Biostrings)

yeast = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL'
human = 'MAAQGEPQVQFKLVLVGDGGTGKTTFVKRHLTGEFEKKYVATLGVEVHPLVFHTNRGPIKFNVWDTAGQEKFGGLRDGYYIQAQCAIIMFDVTSRVTYKNVPNWHRDLVRVCENIPIVLCGNKVDIKDRKVKAKSIVFHRKKNLQYYDISAKSNYNFEKPFLWLARKLIGDPNLEFVAMPALAPPEVVMDPALAAQYEHDLEVAQTTALPDEDDDL'
canine = 'MAAQGEPQVQFKLVLVGDGGTGKTTFVKRHLTGEFEKKYVATLGVEVHPLVFHTNRGPIKFNVWDTAGQEKFGGLRDGYYIQAQCAIIMFDVTSRVTYKNVPNWHRDLVRVCENIPIVLCGNKVDIKDRKVKAKSIVFHRKKNLQYYDISAKSNYNFEKPFLWLARKLIGDPNLEFVAMPALAPPEVVMDPALAAQYEHDLEVAQTTALPDEDDDL'


# compare yeast to human
p <- pairwiseAlignment(pattern = c(yeast), subject = human)
frac_match = nmatch(p)/nchar(human)

#p <- pairwiseAlignment(pattern = c(yeast), subject = canine)
#frac_match = nmatch(p)/nchar(canine)

perc_sequence_identity = pid(p)

# yeast residue numbers are 2 larger than human starting at human 8 / yeast 10
# yeast residue numbers are 2 larger than canine starting at canine 8 / yeast 10 --> same as human!
p@pattern@mismatch@unlistData
p@subject@mismatch@unlistData
p@pattern@mismatch@unlistData - p@subject@mismatch@unlistData