## 03_merge VIRION and Upham phylogeny
## danbeck@ou.edu
## Updated by Briana Betke on 10/29/25

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(ape)
#library(Hmisc) #caps corrected in previous files
#library(plyr)

## load species-level virus data
setwd("/Users/brianabetke/Desktop/bathaus")
virus=read.csv("virus data/viral response_bats VIRION flat.csv")

## load in Upham phylogeny
tree=read.nexus('phylos/MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## load in taxonomy
taxa=read.csv('phylos/taxonomy_mamPhy_5911species.csv',header=T)
taxa=taxa[taxa$ord=="CHIROPTERA",]
taxa$tip=taxa$Species_Name

## trim phylo to bats
tree=keep.tip(tree,taxa$tiplabel)

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## correct VIRION names
virus$species=virus$Host
virus$Host=NULL

# fix genus modifications
virus$species <- gsub("Neoeptesicus|Cnephaeus","Eptesicus",virus$species)

## which names in VIRION aren't in Upham
setdiff(virus$species,taxa$species)

## match manually
virus$species=plyr::revalue(virus$species,
                            c("Artibeus cinereus" = "Dermanura cinereus",
                              "Artibeus glaucus"="Dermanura glaucus",
                              "Artibeus toltecus"="Dermanura toltecus",
                              "Doryrhina cyclops"="Hipposideros cyclops",
                              "Epomophorus dobsonii" = "Epomops dobsonii",
                              "Epomophorus pusillus" = "Micropteropus pusillus",
                              "Glossophaga mutica" = "Glossophaga soricia",
                              "Hypsugo alaschanicus"="Pipistrellus alaschanicus",
                              "Hypsugo pulveratus"="Pipistrellus pulveratus",
                              "Hypsugo savii"="Pipistrellus savii",
                              #"Kerivoula furva" = "Kerivoula titania",
                              "Laephotis capensis"="Neoromicia capensis",
                              "Lyroderma lyra"="Megaderma lyra",
                              "Macronycteris commersonii"="Hipposideros commersoni",
                              "Macronycteris gigas"="Hipposideros gigas",
                              "Macronycteris vittatus"="Hipposideros vittatus",
                              #"Miniopterus africanus"="Miniopterus inflatus",
                              #"Molossus nigricans"="Molossus rufus",
                              "Mops plicatus"="Chaerephon plicatus",
                              "Mops pumilus" = "Chaerephon pumilus",
                              "Murina feae"="Murina aurata",
                              #"Myotis oxygnathus"="Myotis blythii",
                              "Nycticeinops crassulus"="Pipistrellus crassulus",
                              #Otomops harrisoni
                              "Rhinolophus hildebrandti" = "Rhinolophus hildebrandtii",
                              "Perimyotis subflavus"="Pipistrellus subflavus",
                              "Plecotus gaisleri"="Plecotus teneriffae",
                              "Pseudoromicia brunnea"="Neoromicia brunnea",
                              "Pseudoromicia tenuipinnis"="Neoromicia tenuipinnis"
                            ))

## recheck missing
setdiff(virus$species,taxa$species)

## remove missing
virus=virus[!virus$species%in%setdiff(virus$species,taxa$species),]

## do we have duplicate species names after taxonomic matching?
virus[duplicated(virus$species),]

## clean taxonomy file
taxa=taxa[c("species","gen","fam","clade","MSW3_sciName_matched","tip")]

## merge with species
data=merge(taxa,virus,by="species",all.x=T)

## Pseudo absences 
data <- replace_na(data, list(vfam = 0, Virus = 0, dum_virus = 0))

## Zoonotic data
# read in zvirus data
# setwd("/Volumes/BETKE 2021/bathaus/virus data")
zvirus=read_csv("virus data/zoonotic viral response_bats VIRION flat.csv") 

## correct VIRION names
zvirus$species=zvirus$Host
# zvirus$X=NULL
zvirus$Host=NULL

# fix genus modifications
zvirus$species <- gsub("Neoeptesicus|Cnephaeus","Eptesicus", zvirus$species)

## which names in VIRION aren't in Upham
setdiff(zvirus$species,taxa$species)

## match manually
zvirus$species=plyr::revalue(zvirus$species,
                             c("Artibeus cinereus" = "Dermanura cinereus",
                               "Artibeus glaucus"="Dermanura glaucus",
                               "Artibeus toltecus"="Dermanura toltecus",
                               "Doryrhina cyclops"="Hipposideros cyclops",
                               "Epomophorus dobsonii" = "Epomops dobsonii",
                               "Epomophorus pusillus" = "Micropteropus pusillus",
                               "Glossophaga mutica" = "Glossophaga soricia",
                               "Hypsugo alaschanicus"="Pipistrellus alaschanicus",
                               "Hypsugo pulveratus"="Pipistrellus pulveratus",
                               "Hypsugo savii"="Pipistrellus savii",
                               #"Kerivoula furva" = "Kerivoula titania",
                               "Laephotis capensis"="Neoromicia capensis",
                               "Lyroderma lyra"="Megaderma lyra",
                               "Macronycteris commersonii"="Hipposideros commersoni",
                               "Macronycteris gigas"="Hipposideros gigas",
                               "Macronycteris vittatus"="Hipposideros vittatus",
                               #"Miniopterus africanus"="Miniopterus inflatus",
                               #"Molossus nigricans"="Molossus rufus",
                               "Mops plicatus"="Chaerephon plicatus",
                               "Mops pumilus" = "Chaerephon pumilus",
                               "Murina feae"="Murina aurata",
                               #"Myotis oxygnathus"="Myotis blythii",
                               "Nycticeinops crassulus"="Pipistrellus crassulus",
                               #Otomops harrisoni
                               "Rhinolophus hildebrandti" = "Rhinolophus hildebrandtii",
                               "Perimyotis subflavus"="Pipistrellus subflavus",
                               "Plecotus gaisleri"="Plecotus teneriffae",
                               "Pseudoromicia brunnea"="Neoromicia brunnea",
                               "Pseudoromicia tenuipinnis"="Neoromicia tenuipinnis"
                             ))

## recheck missing
setdiff(zvirus$species,taxa$species)

## remove missing
zvirus=zvirus[!zvirus$species%in%setdiff(zvirus$species,taxa$species),]

## do we have duplicate species names after taxonomic matching?
zvirus[duplicated(zvirus$species),]

## merge with species
data=merge(data,zvirus,by="species",all.x=T)

## if na zoonotic viruses, 0 zoonotic viruses
data <- replace_na(data, list(zoo_sp = 0, zfam = 0, dum_zvirus = 0))

## check that there are no instances where zoonotic richness is greater than overall
table(data$zoo_sp > data$Virus)
table(data$zfam > data$vfam)

## write flat file
# setwd("/Volumes/BETKE 2021/bathaus/flat files")
write.csv(data,"flat files/bathaus_virus to phylo backbone.csv")
