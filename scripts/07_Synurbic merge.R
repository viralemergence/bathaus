# 07_ Synurbic merge
# Merging synurbic data to virus COMBINE merged dataset
# babetke@utexas.edu

# Updated October 2025

# clear environment
rm(list=ls())
graphics.off()

# packages
library(tidyverse)

setwd("/Users/brianabetke/Desktop/bathaus")
# taxonomy, citations, virus, synurb skeleton data
virus <- read_csv("flat files/master data_1287 species.csv") %>%
  select(!(...1)) 

# Google sheet data
batroost <- read_csv("Synurbic data/Bat References Spreadsheet.csv") %>%
  select(species, Synurbic, Complete)

# look at how many complete observations you have - should be only yes
table(batroost$Complete)

# merge
merged <- merge(virus, batroost, by = "species", all = TRUE) %>%
  rename(Synurbic = Synurbic.y) %>%
  select(!(Synurbic.x)) %>%
  mutate(fam = ifelse(gen == "Miniopterus", "MINIOPTERIDAE", fam))

# Export
#setwd("/Volumes/BETKE 2021/bathaus/flat files")
write.csv(merged, "flat files/Synurb and filter data Merge.csv", row.names = FALSE)
