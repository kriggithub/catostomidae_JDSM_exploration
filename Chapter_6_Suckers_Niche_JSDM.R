#### K. Riggin May 5, 2025
#### Joint Species Distribution Modeling Traits and Phylogeny

library(ape)


data <- read.csv("suckerspresencedata.csv")
data$elevation <- data$elevation/100

# separate model construction and validation sites
data_C <- data %>% filter(val_cond == "C")
data_V <- data %>% filter(val_cond == "V")



tree <- read.nexus("mito_tree.nex")
plot(tree)

# construct phylogenetic correlation matrix
C <- vcv(tree, model = "Brownian", corr = T)



#influence of traits on suckers




traits <- read.csv("catostomidae_traits.csv")
