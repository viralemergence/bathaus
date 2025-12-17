## 01_virus pull
## pulling all bat viruses from VIRION
## danbeck@ou.edu

# updated to Virion V8 by Briana Betke on 10/24/25

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(vroom)
library(ape)

## load VIRION
# V8 of virion
virion <- vroom(file = "/Users/brianabetke/outputs/15692263/virion.csv.gz")

# update host names to capitalize to be consistent with PREDICT names
virion$Host <- Hmisc::capitalize(virion$Host)

# Pull all chiroptera, including PREDICT which is capitalized
bats <- virion %>% filter(HostOrder == "chiroptera"| HostOrder =="Chiroptera")

# NCBI ratified hosts and viruses
bats <- bats %>% filter(HostNCBIResolved==T & VirusNCBIResolved==T)

## Descriptives of detections
# How many observations missing VirusFamily?
bats %>% drop_na(Host, Virus) %>% count(is.na(VirusFamily)) #298

# Are there more than one detection per pair with missing virus family
bats %>% drop_na(Host, Virus) %>% mutate(ID = paste(Host, Virus)) %>% 
  group_by(ID) %>% tally(is.na(VirusFamily)) %>% filter(n > 0) %>% rename(missing = n) -> miss
# 214 unique pairs with missing virus families

# Pull tally for nonmissing virus families
bats %>% drop_na(Host, Virus) %>% mutate(ID = paste(Host, Virus)) %>% 
  group_by(ID) %>% tally(!is.na(VirusFamily)) %>% rename(non = n) -> non

# build a df to look at how many pairs without virus family also have detections
# with family 
merge(miss, non, by = "ID") -> dat

# total number of detections per pair (not unique detection method)
bats %>% drop_na(Host, Virus) %>% mutate(ID = paste(Host, Virus)) %>% 
  group_by(ID) %>% count() %>% rename(total = n) -> tot

# merge by ID for full comparison
dat <- merge(dat, tot, by = "ID") 
# it looks like none of the pairs have a mixture of detections with virus family and without

# remove
rm(miss,tot,non,dat)

# look at tally of detection methods
bats %>%
  group_by(Host, Virus) %>%
  distinct(DetectionMethod) %>%
  drop_na() %>%
  mutate(dt = 1) %>%
  pivot_wider(names_from = DetectionMethod, values_from = dt) %>%
  replace(is.na(.), 0) %>%
  ungroup() -> det

det %>% 
  mutate(method = rowSums(across(c("Not specified","Isolation/Observation","PCR/Sequencing","Antibodies"))),
         sero = ifelse(Antibodies == 1 & method == 1, 1, 0),
         iso = ifelse(`Isolation/Observation` == 1 & method == 1, 1, 0),
         pcr = ifelse(`PCR/Sequencing` == 1 & method == 1, 1, 0),
         unknown = ifelse(`Not specified` == 1 & method == 1, 1, 0)) -> det

# how many pairs relay solely on each detection method?
det %>% summarise(sero = sum(sero),
                  iso = sum(iso),
                  unknown = sum(unknown),
                  pcr = sum(pcr))

# clean
rm(det)

## Cleaning
# remove no virus family and na
bats <- bats %>% drop_na(Host, Virus) %>% filter(!is.na(VirusFamily))

# Filter to isolation and PCR
bats %>% filter(DetectionMethod == "PCR/Sequencing" | DetectionMethod == "Isolation/Observation") -> bats

# check the detection methods are only PCR or isolation
table(bats$DetectionMethod)

## Calculate family and species richness
aggregate(cbind(VirusFamily,Virus)~Host,data=bats, function(x) length(unique(x))) -> agg

# Rename col
agg %>% rename(vfam = VirusFamily) -> agg

# binary hosting
agg %>% mutate(dum_virus = 1) -> agg

# ## merge to phylo backbone
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
# ## correct VIRION names
# agg$species=agg$Host
# agg$Host=NULL
# 
# ## which names in VIRION aren't in Upham
# setdiff(agg$species,taxa$species)
# 
# # fix genus modifications
# agg$species <- gsub("Neoeptesicus|Cnephaeus","Eptesicus", agg$species)
# 
# # more difs
# setdiff(agg$species,taxa$species)
# 
# # update names
# agg$species=plyr::revalue(agg$species,
#                               c("Artibeus cinereus" = "Dermanura cinereus",
#                                 "Artibeus glaucus"="Dermanura glaucus",
#                                 "Artibeus toltecus"="Dermanura toltecus",
#                                 "Doryrhina cyclops"="Hipposideros cyclops",
#                                 "Epomophorus dobsonii" = "Epomops dobsonii",
#                                 "Epomophorus pusillus" = "Micropteropus pusillus",
#                                 "Glossophaga mutica" = "Glossophaga soricia",
#                                 "Hypsugo alaschanicus"="Pipistrellus alaschanicus",
#                                 "Hypsugo pulveratus"="Pipistrellus pulveratus",
#                                 "Hypsugo savii"="Pipistrellus savii",
#                                 #"Kerivoula furva" = "Kerivoula titania",
#                                 "Laephotis capensis"="Neoromicia capensis",
#                                 "Lyroderma lyra"="Megaderma lyra",
#                                 "Macronycteris commersonii"="Hipposideros commersoni",
#                                 "Macronycteris gigas"="Hipposideros gigas",
#                                 "Macronycteris vittatus"="Hipposideros vittatus",
#                                 #"Miniopterus africanus"="Miniopterus inflatus",
#                                 #"Molossus nigricans"="Molossus rufus",
#                                 "Mops plicatus"="Chaerephon plicatus", 
#                                 "Mops pumilus" = "Chaerephon pumilus",
#                                 "Murina feae"="Murina aurata",
#                                 #"Myotis oxygnathus"="Myotis blythii",
#                                 "Nycticeinops crassulus"="Pipistrellus crassulus",
#                                 #Otomops harrisoni
#                                 "Rhinolophus hildebrandti" = "Rhinolophus hildebrandtii",
#                                 "Perimyotis subflavus"="Pipistrellus subflavus",
#                                 "Plecotus gaisleri"="Plecotus teneriffae",
#                                 "Pseudoromicia brunnea"="Neoromicia brunnea",
#                                 "Pseudoromicia tenuipinnis"="Neoromicia tenuipinnis"
#                               ))
# 
# # check for missing and duplicates
# setdiff(agg$species,taxa$species)
# 
# # clear missing
# agg=agg[!agg$species%in%setdiff(agg$species,taxa$species),]
# 
# # look at duplicates
# agg[duplicated(agg$species),]
# 
# # clean duplicates
# # agg <- agg[!duplicated(fam_agg$species),]
# 
# # merge
# mtax <- merge(taxa, agg, by = "species", all = TRUE)
# 
# # make NAs zero
# mtax <- replace_na(mtax, list(vfam = 0, Virus = 0, dum_virus = 0))

# save 
write.csv(agg,"/Users/brianabetke/Desktop/bathaus/virus data/viral response_bats VIRION flat.csv", row.names = FALSE)
