setwd("~/R/catostomidae_JDSM_exploration/phylogenetic_tree")
mitochondrial <- readxl::read_excel("Table S1. Mitochondrial dataset sample information.xlsx")
# extract accession numbers and match to species names
mitoaccessions <- mitochondrial$...7
accession_species <- mitochondrial$...5

rows <- !is.na(mitoaccessions)
mitoaccessions_clean <- mitoaccessions[rows]
mitoaccessions_clean <- mitoaccessions_clean[-1]
mitoaccessions_clean

accession_species <- accession_species[rows]
accession_species <- accession_species[-1]
accession_species


library(ape)
library(DECIPHER)
library(seqinr)

# pull sequences from GenBank
sequences <- read.GenBank(mitoaccessions_clean)

# convert sequences to fasta format
write.dna(sequences, file = 'mito_sequences', format = 'fasta')

# load in fasta file
fas <- "~/R/catostomidae_JDSM_exploration/phylogenetic_tree/mito_sequences"
seqs <- readDNAStringSet(fas)

# view sequenecs
seqs

# orient and align nucleotides
seqs <- OrientNucleotides(seqs)
aligned <- AlignSeqs(seqs)


# view alignment in browers
BrowseSeqs(aligned, highlight=0)


# write alignment to a new FASTA file
writeXStringSet(aligned, file="mito_aligned.fasta")


# read in aligned data
dna <- read.alignment("mito_aligned.fasta", format = "fasta")

# create a distance matrix for the alignment 
D <- dist.alignment(dna, matrix = "similarity")


# initialize tree
tree <- nj(D)
tree <- ladderize(tree)



# change tree tip labels
tree$tip.label
tree$tip.label <- accession_species
tree$tip.label

# plot tree
plot(tree, cex = 0.6)


# export tree
write.nexus(tree, file = "mito_tree.nex")

