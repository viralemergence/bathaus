# 06_COMBINE merge 
# Merging master google sheet with COMBINE trait dataset
# babetke@utexas.edu

rm(list=ls()) 
graphics.off()

# data wranglin
library(tidyverse)

# Read in the merged data
merge_data <- read_csv("~/Desktop/Bats and Viruses/bathaus/flat files/Synurb and filter data Merge.csv")

# reading chiroptera combine
trait_data <- read_csv("~/Desktop/Bats and Viruses/bathaus/COMBINE Datasets/COMBINE_archives/trait_data_reported.csv") %>%
  filter(order == "Chiroptera") # returns 1342 observations

# write dataset of just Chiroptera
#write.csv(trait_data, "COMBINE Datasets/COMBINE Chiroptera_1342 Species.csv", row.names = FALSE)

# look at differences, make sure that you change names for the combine dataset to the merged
# names in binomial, not in merge
setdiff(trait_data$iucn2020_binomial, merge_data$species)

# might actually be easier to look at both
setdiff(merge_data$species, trait_data$iucn2020_binomial)

# Matching the mismatches
trait_data$iucn2020_binomial <- trait_data$iucn2020_binomial %>% 
  recode("Anoura carishina" = "Anoura canishina",
         "Scotonycteris ophiodon" = "Casinycteris ophiodon",
         "Chiroderma vizzotoi" = "Chiroderma vizottoi",
         "Dermanura azteca" = "Dermanura aztecus",
         "Dermanura cinerea" = "Dermanura cinereus",
         "Dermanura glauca" = "Dermanura glaucus",
         "Dermanura gnoma" = "Dermanura gnomus", 
         "Dermanura rosenbergi" = "Dermanura rosenbergii",
         "Dermanura tolteca" = "Dermanura toltecus",
         "Pteropus leucopterus" = "Desmalopex leucopterus",
         "Diclidurus isabella" = "Diclidurus isabellus",
         "Paremballonura atrata" = "Emballonura atrata",
         "Paremballonura tiavato" = "Emballonura tiavato",
         "Rhyneptesicus nasutus" = "Eptesicus nasutus",
         "Hypsugo affinis" = "Falsistrellus affinis",
         "Macronycteris commersoni" = "Hipposideros commersoni",
         "Macronycteris gigas" = "Hipposideros gigas",
         "Hipposideros pendlebury" = "Hipposideros pendelburyi",
         "Macronycteris thomensis" = "Hipposideros thomensis",
         "Macronycteris vittatus" = "Hipposideros vittatus",
         "Lonchophylla cadenai" = "Hsunycteris cadenai",
         "Lonchophylla pattoni" = "Hsunycteris pattoni",
         "Pipistrellus anchietae" = "Hypsugo anchietae",
         "Lyroderma lyra" = "Megaderma lyra",
         "Gardnerycteris crenulatum" = "Mimon crenulatum", 
         "Gardnerycteris koepckeae"  = "Mimon koepckeae",
         "Miniopterus schreibersii" = "Miniopterus fuliginosus",
         "Miniopterus orianae" = "Miniopterus oceanensis",
         "Ozimops loriae" = "Mormopterus loriae",       
         "Murina lorelieae" = "Murina loreliae",          
         "Lissonycteris angolensis" = "Myonycteris angolensis",
         "Mormopterus kalinowskii" = "Nyctinomops kalinowskii", 
         "Hypsugo alaschanicus" = "Pipistrellus alaschanicus",
         "Hypsugo anthonyi" = "Pipistrellus anthonyi",   
         "Hypsugo arabicus" = "Pipistrellus arabicus",
         "Hypsugo ariel" = "Pipistrellus ariel",       
         "Hypsugo cadornae" = "Pipistrellus cadornae",  
         "Hypsugo eisentrauti" = "Pipistrellus eisentrauti",
         "Hypsugo joffrei" = "Pipistrellus joffrei",     
         "Hypsugo kitcheneri" = "Pipistrellus kitcheneri",
         "Hypsugo lophurus" = "Pipistrellus lophurus",   
         "Hypsugo macrotis" = "Pipistrellus macrotis",
         "Hypsugo musciculus" = "Pipistrellus musciculus",  
         "Hypsugo pulveratus" = "Pipistrellus pulveratus",
         "Hypsugo savii" = "Pipistrellus savii",       
         "Perimyotis subflavus" = "Pipistrellus subflavus",
         "Hypsugo vordermanni" = "Pipistrellus vordermanni",
         "Baeodon alleni" = "Rhogeessa alleni",
         "Baeodon gracilis" = "Rhogeessa gracilis",      
         "Boneia bidens" = "Rousettus bidens",
         "Austronomus australis" = "Tadarida australis",      
         "Chaerephon jobimena" = "Tadarida jobimena",
         "Austronomus kuboriensis" = "Tadarida kuboriensis",          
         "Vampyriscus bidens" = "Vampyressa bidens",
         "Vampyriscus brocki" = "Vampyressa brocki",
         "Vampyriscus nymphaea" = "Vampyressa nymphaea"       
                      ) 
setdiff(merge_data$species, trait_data$iucn2020_binomial)
# went from 89 to 35. 54 bats recoded
setdiff(trait_data$iucn2020_binomial, merge_data$species)
# 138 to 84

# cnames column for synonyms 
merge_data$cnames <- recode(merge_data$species,
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

# create a cnames in trait data for merge
trait_data$cnames <- trait_data$iucn2020_binomial

# merge datasets
trait_merge <- merge(merge_data, trait_data, by = "cnames", all.x = TRUE) %>%
  select(!c(order, family, genus, species.y, phylacine_binomial)) %>%
  rename(species = species.x)

# write after matching
write.csv(trait_merge,"~/Desktop/Bats and Viruses/bathaus/flat files/master data and COMBINE.csv", row.names = FALSE)





