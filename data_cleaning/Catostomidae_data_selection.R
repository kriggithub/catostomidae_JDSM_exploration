#### K. Riggin April 22, 2025
#### This code produces subsetted data for the Catostomidae family as well as trait data
#### Load tables made by Data Prep Code.R prior


####Step 1: Load Packages
library(rfishbase)
library(tidyverse)



####Step 2: Create Trait Dataset
fb_tables()


# ## Get the whole spawning and spawn agg table, joined together:
# spawn <- left_join(fb_tbl("spawning"),  
#                    fb_tbl("spawnagg"), 
#                    relationship = "many-to-many")
spawn <- fb_tbl("spawning")
eco <- fb_tbl("ecosystem")

# Filter taxa down to the desired species
suckers <- load_taxa() %>%  filter(Family == "Catostomidae")


## A "filtering join" (inner join) 
spawn <- spawn %>% inner_join(suckers)
eco <- eco %>% rename(SpecCode = Speccode) %>% inner_join(suckers)


combined <- merge(spawn, suckers)
combined2 <- merge(eco, suckers)









####Step 3: Subset Data sets from Data Prep Code.R to only include Catostomids and Same SITE_ID
colnames(fish_table)
scinames <- common_to_sci(colnames(fish_table))


# Filter scientific names to only include Catostomids 
catscinames <- scinames %>%
  filter(Species %in% suckers$Species)

# Turn into uppercase to match fish_table dataset
catscinames$ComName <- str_to_upper(catscinames$ComName)

# Create subsetted table
fish_table_subset <- fish_table %>%
  select("SITE_ID", any_of(catscinames$ComName))





# Match Environmetal and Reference Table Sites with Fish Table
environmental_table_subset <- environmental_table %>%
  filter(SITE_ID %in% fish_table_subset$SITE_ID)
reference_table_subset <- reference_table %>%
  filter(SITE_ID %in% fish_table_subset$SITE_ID)



# Export Datasets for main analysis
# write.csv(fish_table_subset, file = "fish_count_subset.csv")
# write.csv(environmental_table_subset, file = "environmental_table_subset.csv")
# write.csv(reference_table_subset, file = "reference_table_subset.csv")







