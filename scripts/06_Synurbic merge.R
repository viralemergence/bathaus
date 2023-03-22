# 07_ Synurbic merge
# Merging synurbic data to virus COMBINE merged dataset
# babetke@utexas.edu

# Test script for merging the bat data I have at this point and the virus data. Not COMBINE just yet

# packages
library(tidyverse)

# taxonomy, citations, virus, synurb skeleton data
virus <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/filter data_600 species_editable.csv") %>%
  select(!(...1)) 

# Google sheet data
batroost <- read_csv("/Users/brianabetke/Desktop/Bat References Spreadsheet.xlsx - Sheet1.csv") %>%
  select(species, Synurbic, Complete)

# look at how many complete observations you have - 333
table(batroost$Complete)

# merge
merged <- merge(virus, batroost, by = "species", all = TRUE) %>%
  rename(Synurbic = Synurbic.y) %>%
  select(!(Synurbic.x)) 

# Export
setwd("/Users/brianabetke/Desktop/Bats and Viruses/bathaus")
write.csv(merged, "flat files/Synurb and filter data Merge.csv", row.names = FALSE)
