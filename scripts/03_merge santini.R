## 03_merge santini
## Merging Santini Visitors and Dwellers datasets
## babetke@utexas.edu 

# pachages
library(tidyverse)

# Read in datafiles
dwellers <- read_csv("Santini et al 2018/Chiroptera_dwellers.csv", col_names = TRUE)
visitors <- read_csv("Santini et al 2018/Chiroptera_visitors.csv", col_names = TRUE)

# Adding column to dwellers called UrbanStatus to indicate if dweller or neither
d <- dwellers %>% mutate(UrbanStatus = if_else(Synurbic <= 0, "neither", "dweller"))

# Doing the same for visitors
v <- visitors %>% mutate(UrbanStatus = if_else(Synurbic <= 0, "neither", "visitor"))

# Combine the data by rbind
u <- unique(rbind(d,v))

# change value of Lasiurus borealis from 1 to 0
u[53,10] <- 0
View(u) # check that the value was replaced

# View table of UrbanStatus variable to confirm number of dwellers (11) and visitors (6)
table(u$UrbanStatus)

# Write as csv
write.csv(u, "Santini et al 2018/Santini_dwellers and visitors.csv")


