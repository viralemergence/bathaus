# Zoonotic Virus to merged PanTHERIA data for full dataset
# babetke@utexas.edu

# Clear enviornment
rm(list = ls())

# libraries
library(tidyverse)
library(Hmisc) # capitalize()

# Read in PanTHERIA merge data
setwd("/Users/brianabetke/Desktop/Bats and Viruses/bathaus")
data <- read_csv("flat files/master data and PanTHERIA.csv")

# read in zvirus data
zvirus <- read_csv("virus data/zoonotic viral response_bats VIRION flat.csv") %>% 
  select(-...1)

# need to capitalize genus
zvirus <- zvirus %>% 
  mutate(species = capitalize(Host)) %>%
  select(-Host)

# Look at differences 
setdiff(zvirus$species, data$species)

# match - matches from the overall virus pull
zvirus$species <- recode(zvirus$species,
                          "Aeorestes cinereus"="Lasiurus cinereus",
                          "Aeorestes egregius"="Lasiurus egregius",
                          "Afronycteris nana"="Neoromicia nana",
                          "Antrozous dubiaquercus"="Bauerus dubiaquercus",
                          "Artibeus cinereus"="Dermanura cinereus",
                          "Artibeus glaucus"="Dermanura glaucus",
                          "Artibeus phaeotis"="Dermanura phaeotis",
                          "Artibeus toltecus"="Dermanura toltecus",
                          "Chaerephon leucogaster"="Chaerephon pumilus",
                          "Chaerephon pusillus"="Chaerephon pumilus",
                          "Dasypterus ega"="Lasiurus ega",
                          "Dasypterus intermedius"="Lasiurus intermedius",
                          "Dasypterus xanthinus"="Lasiurus xanthinus",
                          "Dobsonia andersoni"="Dobsonia anderseni",
                          "Dobsonia magna"="Dobsonia moluccensis",
                          "Doryrhina cyclops"="Hipposideros cyclops",
                          "Eptesicus regulus"="Vespadelus regulus",
                          "Eptesicus vulturnus"="Vespadelus vulturnus",
                          "Gardnerycteris crenulatum"="Mimon crenulatum",
                          "Hipposideros cf. ruber"="Hipposideros ruber",
                          "Hipposideros terasensis"="Hipposideros armiger",
                          "Hypsugo alaschanicus"="Pipistrellus alaschanicus",
                          "Hypsugo pulveratus"="Pipistrellus pulveratus",
                          "Hypsugo savii"="Pipistrellus savii",
                          "Laephotis capensis"="Neoromicia capensis",
                          "Lissonycteris angolensis"="Myonycteris angolensis",
                          "Macronycteris commersoni"="Hipposideros commersoni",
                          "Macronycteris gigas"="Hipposideros gigas",
                          "Macronycteris vittata"="Hipposideros vittatus",
                          "Miniopterus africanus"="Miniopterus inflatus",
                          "Miniopterus orianae"="Miniopterus oceanensis",
                          "Molossus ater"="Molossus rufus",
                          "Myotis myotis/blythii"="Myotis myotis",
                          "Myotis oxygnathus"="Myotis blythii",
                          "Myotis ricketti"="Myotis pilosus",
                          "Neoromicia brunneus"="Neoromicia brunnea",
                          "Neoromicia somalicus" = "Neoromicia somalica", # added in zvirus
                          "Nyctalus velutinus"="Nyctalus plancyi",
                          "Parahypsugo crassulus"="Pipistrellus crassulus",
                          "Perimyotis subflavus"="Pipistrellus subflavus",
                          "Plecotus gaisleri"="Plecotus teneriffae",
                          "Pteronotus alitonus"="Pteronotus parnellii",
                          "Pteronotus rubiginosus"="Pteronotus parnellii",
                          "Rhinolophus blythi"="Rhinolophus lepidus",
                          #"Rhinolophus cornutus"="",
                          "Rhinolophus hildebrandti"="Rhinolophus hildebrandtii",
                          "Rhinolophus lobatus"="Rhinolophus landeri",
                          "Rhinolophus monoceros"="Rhinolophus pusillus",
                          "Rhinolophus rhodesiae"="Rhinolophus simulator")


# looks like there are some duplicates
zvirus[duplicated(zvirus$species), ]

# aggregate the duplicates
z <- zvirus %>% 
  group_by(species) %>%
  dplyr::summarize(zvirus = sum(zvirus)) # conflict with Hmis

# check again for dulplicates
z[duplicated(z$species), ]
zvirus <- z
rm(z)

# merge 
data <- merge(data, zvirus, by = "species", all.x = T)

# anything NA to be pseudoabsence 
data <- data %>% 
  mutate(zvirus = replace_na(zvirus, 0))

# check that there are no NAs
sum(is.na(data$zvirus))
sum(!is.na(data$zvirus))

# save csv
write.csv(data, "flat files/zoonotic virus to merged PanTHERIA.csv", row.names = FALSE)


