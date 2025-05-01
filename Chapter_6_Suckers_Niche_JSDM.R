#### K. Riggin April 22, 2025
#### Single Species Distribution Modeling
#### Chosen species: White Sucker (present in 199 sites with max count of 169)

library(ape)


tree <- read.nexus("mito_tree.nex")
plot(tree)

# construct phylogenetic correlation matrix
C <- vcv(tree, model = "Brownian", corr = T)



traits <- read.csv("catostomidae_traits.csv")
