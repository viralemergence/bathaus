# 08_ PanTHERIA merge
# Merging the geographic variables from PanTHERIA to master and COMBINE
# babetke@utexas.edu

# Clear environment
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)

# Read in the merged dataset
setwd("/Users/brianabetke/Desktop/Bats and Viruses/bathaus")
data <- read_csv("flat files/master data and COMBINE.csv") %>% 
  select(-cnames)

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
                            "Artibeus watsoni" = "Dermanura watsoni"
                            ) 
# Check again 
setdiff(data$species, GR$MSW05_Binomial)

# Cnames if I need it for synonyms
data$cnames <- recode(data$species,
                      "Chaerephon pumilus" = "Chaerephon leucogaster",
                      "Chaerephon chapini" = "Chaerephon shortridgei",
                      "Neoromicia matroka" = "Eptesicus matroka",
                      "Carollia brevicauda" = "Carollia colombiana",
                      "Pipistrellus ariel" = "Hypsugo bodenheimeri",
                      "Murina ussuriensis" = "Murina silvatica",
                      "Micronycteris minuta" = "Micronycteris homezi",
                      "Myotis ikonnikovi" = "Myotis yesoensis",
                      "Myotis ikonnikovi" = "Myotis hosonoi",
                      "Myotis blythii" = "Myotis oxygnathus",
                      "Myotis ikonnikovi" = "Myotis ozensis",
                      "Pteropus alecto" = "Pteropus banakrisi",
                      "Dermanura phaeotis" = "Artibeus phaeotis"
                      )

# check to see if they all appear properly 
c <- data %>% filter(species != cnames) %>% select(species,cnames)
rm(c)

# create a cnames in trait data for merge
GR$cnames <- GR$MSW05_Binomial

# merge datasets
data_GR <- merge(data, GR, by = "cnames", all.x = TRUE) %>%
  select(!c(cnames,MSW05_Binomial))

# which one is duplicated? 
data_GR[duplicated(data_GR$cnames),]
# Dermanura phaeotis - should actually be in cnames

# write after matching
write.csv(data_GR,"flat files/master data and PanTHERIA.csv", row.names = FALSE)

# values for the myotis
myotis <- PanTHERIA %>% 
  filter(MSW05_Binomial %in% c("Myotis yesoensis", "Myotis hosonoi", "Myotis ozensis"))
