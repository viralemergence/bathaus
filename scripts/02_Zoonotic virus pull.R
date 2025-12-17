## 02_Zoonotic Virus Pull
## pulling all zoonotic viruses from VIRION
## danbeck@ou.edu
## Updated 10/27/25 by Briana Betke (update virion to V8)

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(vroom)
library(ape)

## load VIRION
#setwd("~/Desktop/Virion")
virion <- vroom(file = "/Users/brianabetke/outputs/15692263/virion.csv.gz")

# update host names to capitalize genus in PREDICT names
virion$Host <- Hmisc::capitalize(virion$Host)

# Pull all chiroptera, including PREDICT which is capitalized
virion <- virion %>% filter(HostOrder == "chiroptera"| HostOrder =="Chiroptera"| Host == "Homo sapiens")

#NCBI ratified hosts and viruses
virion <- virion %>% filter(HostNCBIResolved==T & VirusNCBIResolved==T)

# Drop NAs
virion <- virion %>% drop_na(Host, Virus) %>% filter(!is.na(VirusFamily))

# trim to PCR and Isolation only
virion <- virion %>% filter(DetectionMethod == "PCR/Sequencing" | DetectionMethod == "Isolation/Observation")

# trim out duplicate associations
virion <- virion %>% mutate(unique = paste(Host,Virus)) %>% filter(!duplicated(unique))

# double check
length(unique(virion$unique))

# filter by humans and bats
humans <- virion %>% filter(Host == "Homo sapiens")
bats <- virion %>% filter(HostOrder == "chiroptera"| HostOrder =="Chiroptera")

# pull out human viruses
hvir <- humans$Virus

# calculate species
bats <- bats %>% mutate(zoo_sp = ifelse(Virus %in% hvir,1,0))

# aggregate 1s
sp_agg <- aggregate(zoo_sp~Host, data = bats, function(x) sum(x))

## Calculate richness values
hb <- bats %>% filter(Virus %in% hvir)

## calculate family richness
zfam_agg <- aggregate(VirusFamily~Host, data=hb, function(x) length(unique(x)))

# rename
zfam_agg <- zfam_agg %>% rename(zfam = VirusFamily)

# aggreate species and family
agg <- merge(sp_agg, zfam_agg, by ="Host", all = TRUE)

# clean out everything but zoonotic
rm(list = ls()[!ls() %in% c("agg")])

# add binary hosting
agg <- agg %>% mutate(dum_zvirus = ifelse(zoo_sp >= 1,1,0)) 

# ### merge to phylo backbone
# # load upham
# tree=read.nexus(here::here("phylos",'MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre'))
# 
# ## load in taxonomy
# taxa=read.csv(here::here("phylos",'taxonomy_mamPhy_5911species.csv'), header=T) 
# taxa=taxa[taxa$ord=="CHIROPTERA",]
# taxa$tip=taxa$Species_Name
# 
# ## trim phylo to bats
# tree=keep.tip(tree,taxa$tiplabel)
# 
# ## fix tip
# tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
# taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))
# 
# ## capitalize virion host names
# agg$species=Hmisc::capitalize(agg$Host)
# agg$Host=NULL
# 
# ## which names in VIRION aren't in Upham
# setdiff(agg$species,taxa$species)
# 
# # fix genus modifications
# agg$species <- gsub("Neoeptesicus|Cnephaeus","Eptesicus", agg$species)
# 
# # more difs
# setdiff(agg$species,agg$species)
# 
# # update names
# agg$species=plyr::revalue(agg$species,
#                                c("Artibeus cinereus" = "Dermanura cinereus",
#                                  "Artibeus glaucus"="Dermanura glaucus",
#                                  "Artibeus toltecus"="Dermanura toltecus",
#                                  "Doryrhina cyclops"="Hipposideros cyclops",
#                                  "Epomophorus dobsonii" = "Epomops dobsonii",
#                                  "Epomophorus pusillus" = "Micropteropus pusillus",
#                                  "Glossophaga mutica" = "Glossophaga soricia",
#                                  "Hypsugo alaschanicus"="Pipistrellus alaschanicus",
#                                  "Hypsugo pulveratus"="Pipistrellus pulveratus",
#                                  "Hypsugo savii"="Pipistrellus savii",
#                                  #"Kerivoula furva" = "Kerivoula titania",
#                                  "Laephotis capensis"="Neoromicia capensis",
#                                  "Lyroderma lyra"="Megaderma lyra",
#                                  "Macronycteris commersonii"="Hipposideros commersoni",
#                                  "Macronycteris gigas"="Hipposideros gigas",
#                                  "Macronycteris vittatus"="Hipposideros vittatus",
#                                  #"Miniopterus africanus"="Miniopterus inflatus",
#                                  #"Molossus nigricans"="Molossus rufus",
#                                  "Mops plicatus"="Chaerephon plicatus", 
#                                  "Mops pumilus" = "Chaerephon pumilus",
#                                  "Murina feae"="Murina aurata",
#                                  #"Myotis oxygnathus"="Myotis blythii",
#                                  "Nycticeinops crassulus"="Pipistrellus crassulus",
#                                  #Otomops harrisoni
#                                  "Rhinolophus hildebrandti" = "Rhinolophus hildebrandtii",
#                                  "Perimyotis subflavus"="Pipistrellus subflavus",
#                                  "Plecotus gaisleri"="Plecotus teneriffae",
#                                  "Pseudoromicia brunnea"="Neoromicia brunnea",
#                                  "Pseudoromicia tenuipinnis"="Neoromicia tenuipinnis"
#                                ))
# 
# # check for missing and duplicates
# setdiff(agg$species,taxa$species)
# 
# # clear missing
# agg <- agg[!agg$species%in%setdiff(agg$species,taxa$species),]
# 
# # look at duplicates
# agg[duplicated(agg$species),]
# # filter(zoonotic, duplicated(species))
# 
# # # clean duplicates
# # zoonotic <- zoonotic[!duplicated(zoonotic$species),]
# 
# # merge
# mtax <- merge(taxa, agg, by = "species", all = TRUE)
# 
# # make NAs zero
# mtax <- replace_na(mtax, list(zoo_sp = 0, zfam = 0))
# 
# # add binary hosting
# mtax <- mtax %>% mutate(dum_zvirus = ifelse(zoo_sp >= 1, 1, 0))

# save 
write.csv(agg,"/Users/brianabetke/Desktop/bathaus/virus data/zoonotic viral response_bats VIRION flat.csv", row.names = FALSE)
