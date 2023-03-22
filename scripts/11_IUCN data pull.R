# 11_IUCN data pull and merge.R
# Data pull for IUCN data by iucn and merge

# Script for data cleaning
rm(list=ls()) 
graphics.off()

# packages
library(tidyverse)
library(rredlist)

# Update synurbic data
## load in previous trait data
zoonotic <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/zoonotic virus to merged PanTHERIA.csv")

# merge updated synurbic data if needed
update <- read_csv("~/Desktop/Synurbic_Bats/synurbat/raw/Bat References Spreadsheet.csv") %>%
  select(species, Synurbic, Complete)# pull from google drive

# check for complete - finished data should be only Yes
table(update$Complete)

# names of bats that are left to revisit - should return nothing when revisits are done
update[update$Complete == "revisit",]$species

# merge updated synurbic data
merged <- merge(zoonotic, update, by = "species", all = TRUE)

merged %>%
  rename(Synurbic = Synurbic.y,
         Complete = Complete.y) %>%
  select(!c(Synurbic.x, Complete.x)) -> merged

# remove NAs in iucn2020 to get only combine names, then do the same matching you did for COMBINE
cdata <- merged %>% drop_na(iucn2020_binomial)

## rredlist scrape - key
bapi="b0f54860abaf1e70f96157345f999eda91e1652e9f199c9d83712c5f64ea6dde"

## run IUCN function
rdata=lapply(cdata$iucn2020_binomial,function(x){
  res=rl_search(name=x,key=bapi)$result; 
  print(x)
  return(res)
}
)

## rbind
rdata=do.call(rbind.data.frame,rdata)
rdata$iucn2020_binomial=rdata$scientific_name

# duplicates
rdata$iucn2020_binomial[duplicated(rdata$iucn2020_binomial)]
# [1] "Dermanura watsoni"      "Harpiocephalus harpia"  "Lonchophylla thomasi"  
# [4] "Lophostoma carrikeri"   "Molossus coibensis"     "Myotis formosus"       
# [7] "Myotis simus"           "Natalus mexicanus"      "Natalus stramineus"    
# [10] "Nyctophilus corbeni"    "Pipistrellus kuhlii"    "Pteropus chrysoproctus"
# [13] "Pteropus pelewensis"    "Rhinolophus borneensis" "Triaenops rufus" 

# remove them
rdata <- unique(rdata)

# cnames (recode names in our dataset)
merged$cnames <- recode(merged$species,
                            "Dermanura incomitatus"="Dermanura watsoni",
                            "Harpiocephalus mordax"="Harpiocephalus harpia",
                            "Hsunycteris thomasi" = "Lonchophylla thomasi",
                            "Lophostoma aequatorialis" = "Lophostoma occidentalis",
                            "Lophostoma yasuni" = "Lophostoma carrikeri",
                            "Miniopterus fuliginosus" ="Miniopterus schreibersii",
                            "Molossus barnesi" = "Molossus coibensis",
                            "Myotis flavus" = "Myotis formosus",
                            "Myotis midastactus" = "Myotis simus",
                            "Natalus lanatus" = "Natalus mexicanus",
                            "Natalus saturatus" = "Natalus stramineus",
                            "Nyctophilus timoriensis" = "Nyctophilus corbeni",
                            "Paracoelops megalotis" = "Hipposideros Pomona",
                            "Pipistrellus deserti" = "Pipistrellus kuhlii",
                            "Pteropus argentatus" = "Pteropus chrysoproctus",
                            "Pteropus yapensis" = "Pteropus pelewensis",
                            "Triaenops menamena" = "Triaenops rufus",
                            "Rhinolophus chaseni" = "Rhinolophus borneensis"
)


# create cnames in cdata
rdata$cnames <- rdata$scientific_name

# merge (should be 1287 total obs)
cname_merge <- merge(merged, rdata[c("cnames","iucn2020_binomial","category","population_trend")], by = "cnames", all.x = TRUE) %>%
  select(!iucn2020_binomial.y) %>%
  rename(iucn2020_binomial = iucn2020_binomial.x)

## missing
table(is.na(cname_merge$category), cname_merge$Synurbic)

# write as csv for merge file
write.csv(cname_merge, "/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/IUCN data merge to zoonotic.csv")
