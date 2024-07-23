# 09_PanTHERIA merge
# Merging the geographic variables from PanTHERIA to master and COMBINE
# babetke@utexas.edu

# Clear environment
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)

# Read in the merged dataset
setwd("/Volumes/BETKE 2021/bathaus")
data <- read_csv("flat files/master data and COMBINE.csv") #%>% 
  #select(-cnames)

# reading in PanTHERIA
PanTHERIA <- read.delim("traits/PanTHERIA_1-0_WR05_Aug2008.txt")

# Grab the geographic variables and species names from PanTHERIA to merge
GR <- PanTHERIA %>%
  filter(MSW05_Order == "Chiroptera") %>%
  select(MSW05_Binomial, X26.1_GR_Area_km2:X30.2_PET_Mean_mm) %>%
  na_if(-999) # change -999 to NA

# look at differences with setdiff()
# Everything in virus, not in PanTHERIA
setdiff(data$species, GR$MSW05_Binomial)
# Everything in PanTHERIA that isn't in virus
setdiff(GR$MSW05_Binomial, data$species)

# # values for 3 myotis synonyms
# myotis <- PanTHERIA %>% 
#   filter(MSW05_Binomial %in% c("Myotis yesoensis", "Myotis hosonoi", "Myotis ozensis"))

# Matching the mismatches
GR$MSW05_Binomial <- recode(GR$MSW05_Binomial, 
                            "Hypsugo alaschanicus" = "Pipistrellus alaschanicus",
                            "Hypsugo anthonyi" = "Pipistrellus anthonyi",
                            "Hypsugo arabicus" = "Pipistrellus arabicus",
                            "Hypsugo ariel" = "Pipistrellus ariel",
                            "Hypsugo cadornae" = "Pipistrellus cadornae",
                            "Hypsugo crassulus" = "Pipistrellus crassulus",
                            "Hypsugo eisentrauti" = "Pipistrellus eisentrauti",
                            "Hypsugo imbricatus" = "Pipistrellus imbricatus",
                            "Hypsugo joffrei" = "Pipistrellus joffrei",
                            "Hypsugo kitcheneri" = "Pipistrellus kitcheneri",
                            "Hypsugo lophurus" = "Pipistrellus lophurus",
                            "Hypsugo macrotis" = "Pipistrellus macrotis",
                            "Hypsugo musciculus" = "Pipistrellus musciculus",
                            "Hypsugo pulveratus" = "Pipistrellus pulveratus",
                            "Hypsugo savii" = "Pipistrellus savii",
                            "Hypsugo vordermanni" = "Pipistrellus vordermanni",
                            "Lissonycteris angolensis" = "Myonycteris angolensis",
                            "Murina grisea" = "Harpiola grisea",
                            "Myotis abei" = "Myotis petax",
                            "Myotis ricketti" = "Myotis pilosus",
                            "Neoromicia brunneus" = "Neoromicia brunnea",
                            "Neoromicia nanus" = "Neoromicia nana",
                            "Neoromicia somalicus" = "Neoromicia somalica",
                            "Pteralopex acrodonta" = "Mirimiri acrodonta",
                            "Pteropus insularis" = "Pteropus pelagicus",
                            "Pipistrellus hesperus" = "Parastrellus hesperus",
                            "Plecotus alpinus" = "Plecotus macrobullaris",
                            "Pteropus leucopterus" = "Desmalopex leucopterus",
                            "Scotonycteris ophiodon" = "Casinycteris ophiodon", 
                            "Sturnira thomasi" = "Sturnira angeli",
                            "Triaenops auritus" = "Paratriaenops auritus",
                            "Triaenops furculus" = "Paratriaenops furculus",
                            "Artibeus anderseni" = "Dermanura anderseni",
                            "Artibeus aztecus" = "Dermanura aztecus",
                            "Artibeus cinereus" = "Dermanura cinereus",
                            "Artibeus glaucus" = "Dermanura glaucus",
                            "Artibeus gnomus"= "Dermanura gnomus",
                            "Artibeus incomitatus" = "Dermanura incomitatus",
                            "Artibeus toltecus" = "Dermanura toltecus",
                            "Artibeus watsoni" = "Dermanura watsoni",
                            "Eptesicus matroka" = "Neoromicia matroka",
                            "Artibeus phaeotis" = "Dermanura phaeotis"
                            ) 
# Check again 
setdiff(data$species, GR$MSW05_Binomial)
# went from 237 to 193

# assign traits to synonyms
data$pcnames <- recode(data$species,
                             "Carollia sowelli" = "Carollia brevicauda",
                             "Dermanura incomitatus" = "Dermanura watsoni",
                             "Harpiocephalus mordax" = "Harpiocephalus harpia",
                             "Hsunycteris thomasi" = "Lonchophylla thomasi",
                             "Lophostoma aequatorialis" = "Lophostoma occidentalis",
                             "Lophostoma yasuni" = "Lophostoma carrikeri",
                             "Miniopterus fuliginosus" ="Miniopterus schreibersii",
                             "Molossus barnesi" = "Molossus coibensis",
                             "Myotis flavus" = "Myotis formosus",
                             "Myotis midastactus" = "Myotis simus",
                             "Natalus mexicanus" = "Natalus stramineus",
                             "Natalus lanatus" = "Natalus stramineus",
                             "Natalus saturatus" = "Natalus stramineus",
                             "Nyctophilus corbeni" = "Nyctophilus timoriensis",
                             "Paracoelops megalotis" = "Hipposideros Pomona",
                             "Pipistrellus deserti" = "Pipistrellus kuhlii",
                             "Pteropus argentatus" = "Pteropus chrysoproctus",
                             "Pteropus yapensis" = "Pteropus pelewensis",
                             "Triaenops menamena" = "Triaenops rufus",
                             "Rhinolophus chaseni" = "Rhinolophus borneensis",
                             "Triaenops rufus" = "Triaenops persicus",
                             "Myotis aelleni" = "Myotis chiloensis",
                             "Myotis aurascens" = "Myotis hajastanicus"
)

# check to see if they all appear properly 
c <- data %>% filter(species != pcnames) %>% select(species,pcnames)
rm(c)

# create a cnames in trait data for merge
GR$pcnames <- GR$MSW05_Binomial

# merge datasets, don't remove pcnames
data_GR <- merge(data, GR, by = "pcnames", all.x = TRUE) %>%
  select(!c(MSW05_Binomial))

# check duplicates for pcnames
d <- data_GR[duplicated(data_GR$pcnames),]
rm(d)

# write after matching
write.csv(data_GR,"flat files/master data and PanTHERIA.csv", row.names = FALSE)