## 05_Taxonomy Merge
## Merging combined Santini with Taxonomy data

library(tidyverse)

# read in csv
santini <- read_csv("Santini et al 2018/Santini_dwellers and visitors.csv", col_names = TRUE) %>%
  select(-X1) %>% #remove X1 column
  rename(Species_Name = Species_ph) # change name to match taxonomy file for merging

# Updating Neoromicia_nanus to Neoromicia_nana
santini[48,3] <- "Neoromicia_nana"
santini[48,3]

tax_bats <- read_csv("phylos/taxonomy_mamPhy_5911species.csv", col_names = TRUE) %>%
  filter(ord == "CHIROPTERA") # filter only for bats

# Merge by species name?
f <- merge(tax_bats, santini, on = "Species_Name", all=T)
write_csv(f, "phylos/Santini and Taxonomy merge.csv")





