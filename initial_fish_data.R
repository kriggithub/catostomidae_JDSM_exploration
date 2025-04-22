library(rfishbase)
library(tidyverse)


fb_tables()


# ## Get the whole spawning and spawn agg table, joined together:
# spawn <- left_join(fb_tbl("spawning"),  
#                    fb_tbl("spawnagg"), 
#                    relationship = "many-to-many")
fishdata <- read.csv("")


spawn <- fb_tbl("spawning")
eco <- fb_tbl("ecosystem")

# Filter taxa down to the desired species
suckers <- load_taxa() %>%  filter(Family == "Catostomidae")


## A "filtering join" (inner join) 
spawn <- spawn %>% inner_join(suckers)
eco <- eco %>% rename(SpecCode = Speccode) %>% inner_join(suckers)


combined <- merge(spawn, suckers)
combined2 <- merge(eco, suckers)

common_to_sci()




# build phylogeny and read it in
library(ape)


