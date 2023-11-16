# 07_ Synurbic merge
# Merging synurbic data to virus COMBINE merged dataset
# babetke@utexas.edu

# updated virus data 11/14/23

# clear environment
rm(list=ls())
graphics.off()

# packages
library(tidyverse)

# taxonomy, citations, virus, synurb skeleton data
virus <- read_csv("/Volumes/BETKE 2021/bathaus/flat files/master data_1287 species.csv") %>%
  select(!(...1)) 

# Google sheet data
batroost <- read_csv("/Volumes/BETKE 2021/bathaus/Synurbic data/Bat References Spreadsheet.csv") %>%
  select(species, Synurbic, Complete)

# look at how many complete observations you have - 333
table(batroost$Complete)

# merge
merged <- merge(virus, batroost, by = "species", all = TRUE) %>%
  rename(Synurbic = Synurbic.y) %>%
  select(!(Synurbic.x)) %>%
  mutate(fam = ifelse(gen == "Miniopterus", "MINIOPTERIDAE", fam))

# Export
setwd("/Volumes/BETKE 2021/bathaus/flat files")
write.csv(merged, "Synurb and filter data Merge.csv", row.names = FALSE)
