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
# Version 8 (https://zenodo.org/records/15733485)
virion <- vroom(file = "/Users/brianabetke/Desktop/bathaus/Verena data/virion v8.csv.gz")

# update host names to capitalize genus in PREDICT names
virion$Host <- Hmisc::capitalize(virion$Host)

# Pull all chiroptera and human observations
virion <- virion %>% filter(HostOrder == "chiroptera"| HostOrder =="Chiroptera" | Host == "Homo sapiens")

#NCBI ratified hosts and viruses
virion <- virion %>% filter(HostNCBIResolved==T & VirusNCBIResolved==T)

# Drop NAs
virion <- virion %>% drop_na(Host, Virus)

# trim to PCR and Isolation only
virion <- virion %>% filter(DetectionMethod == "PCR/Sequencing" | DetectionMethod == "Isolation/Observation")

# trim out duplicate associations
virion <- virion %>% mutate(unique = paste(Host,Virus)) %>% filter(!duplicated(unique))

# double check
length(unique(virion$unique))

# clean out names ending in sp.
virion %>% filter(!grepl("sp\\.", Host)) -> virion
virion %>% filter(!grepl("sp\\.", Virus)) -> virion 

# pull and clean bats
bats <- virion %>% filter(HostOrder == "chiroptera"| HostOrder =="Chiroptera")

# check that there is the same number of nonratified names
bats %>% filter(ICTVRatified == F) %>% pull(Virus) %>% unique() %>% length()

bats %>% filter(ICTVRatified == F) %>% pull(Virus) %>% unique() -> nonrat

# read in the renames
renames <- read_csv("/Users/brianabetke/Desktop/bathaus/virus data/Virus reconciliation - renaming.csv")

# check for duplicates
table(duplicated(renames$`virion names`))
# FALSE 
# 1173 

# table
table(renames$status)
# abolished not in ictv     renamed 
#         1        1129          43 

# matches - all 1173
dplyr::intersect(nonrat, renames$`virion names`) %>% length()

# look at the renames
renames %>% filter(status == "renamed") %>% select(`virion names`,`ictv names`) -> renames_df

# make cases consistent with old names
renames_df$`ictv names` <- tolower(renames_df$`ictv names`)

# recode virus names 
bats$Virus <- recode(bats$Virus, 
                     "adenovirus desmodus35011" = "mastadenovirus desmodi",
                     "alphacoronavirus 1" = "alphacoronavirus suis",
                     "avian coronavirus" = "gammacoronavirus galli",
                     "bat associated rhabdovirus 1" = "betaplatrhavirus pipistrellus",
                     "bat associated rhabdovirus 2" = "betaplatrhavirus abramus",
                     "bat associated rhabdovirus 3" = "betaplatrhavirus robustula",
                     "bat coronavirus hku10" = "alphacoronavirus rousetti",
                     "bat cyclovirus" = "cyclovirus flagermus",
                     "bat hepatitis virus" = "orthohepadnavirus longidigiti",
                     "bat ledantevirus 2" = "ledantevirus hipposideros",
                     "bat mastadenovirus" = "mastadenovirus asiense",
                     "bat mastadenovirus a" = "mastadenovirus musauriti",
                     "bat parvovirus" = "champhamavirus chiropteran3",
                     "betacoronavirus 1" = "betacoronavirus gravedinis",
                     "bombali ebolavirus" = "orthoebolavirus bombaliense",
                     "broome orthoreovirus" = "orthoreovirus broomense", 
                     "eastern equine encephalitis virus" = "alphavirus eastern",
                     "eptesicus fuscus gammaherpesvirus" = "patagivirus vespertilionidgamma3",
                     "hepatitis b virus" = "orthohepadnavirus hominoidei",
                     "human coronavirus 229e" = "alphacoronavirus chicagoense",
                     "jeilongvirus madagascarense" = "parajeilongvirus madagascarense",
                     "mamastrovirus 1" = "mamastrovirus hominis",
                     "mamastrovirus 12" = "mamastrovirus vespertilionis",
                     "mamastrovirus 19" = "mamastrovirus guangxiense",
                     "mamastrovirus 9" = "mamastrovirus virginiaens",
                     "mammalian orthoreovirus" = "orthoreovirus mammalis",
                     "marburg marburgvirus" = "orthomarburgvirus marburgense",
                     "middle east respiratory syndrome-related coronavirus" = "betacoronavirus cameli",
                     "miniopterus bat coronavirus hku8" = "alphacoronavirus miniopteri",
                     "mischivirus a" = "mischivirus acombewi",
                     "murine leukemia virus" = "gammaretrovirus murleu",
                     "nipah henipavirus" = "henipavirus nipahense",
                     "orthohantavirus robinaense" = "mobatvirus robinaense",
                     "pipistrellus bat coronavirus hku5" = "betacoronavirus pipistrelli",
                     "pipistrellus kuhlii coronavirus 3398" = "alphacoronavirus pipistrelli",
                     "porcine epidemic diarrhea virus" = "alphacoronavirus porci",
                     "rhinolophus bat coronavirus hku2" = "alphacoronavirus rhinolophi",
                     "rousettus bat coronavirus hku9" = "betacoronavirus rousetti",
                     "scotophilus bat coronavirus 512" = "alphacoronavirus scotophili",
                     "severe acute respiratory syndrome-related coronavirus" = "betacoronavirus pandemicum",
                     "tylonycteris bat coronavirus hku4" = "cetacoronavirus tylonycteridis",
                     "venezuelan equine encephalitis virus" = "alphavirus venezuelan",
                     "western equine encephalitis virus" = "alphavirus western"
)

# update the name for pipestrellus kuhlii parvovirus
bats <- bats %>% mutate(Virus = ifelse(Virus == "champhamavirus chiropteran3" & Host == "Pipistrellus kuhlii", 
                                       "protoparvovirus chiropteran2", Virus))

# change the newly ratified names to to ICTVRatified == TRUE
new_names <- renames_df$`ictv names` 
bats <- bats %>% mutate(ICTVRatified = ifelse(Virus %in% new_names, TRUE, ICTVRatified)) 

# removing the species that no matches were found in genbank for alphacoronavirus rhinolophi
bats <- bats %>% mutate(ICTVRatified = ifelse(Virus == "alphacoronavirus rhinolophi" & Host == "Aselliscus stoliczkanus", 
                                              FALSE, ICTVRatified))

# how many could we gain with just these updates? 202 observations, 57 more from the first round of look ups
table(bats$ICTVRatified)

# filter out remanining non ratified
bats %>% filter(ICTVRatified == T) -> bats

# now, lets pull human viruses
humans <- virion %>% filter(Host == "Homo sapiens")

# see if any of the old names match 
dplyr::intersect(renames_df$`virion names`, humans$Virus)

# update those
# recode virus names 
humans$Virus <- recode(humans$Virus, 
                     "alphacoronavirus 1" = "alphacoronavirus suis",
                     "avian coronavirus" = "gammacoronavirus galli",
                     "betacoronavirus 1" = "betacoronavirus gravedinis",
                     "eastern equine encephalitis virus" = "alphavirus eastern",
                     "hepatitis b virus" = "orthohepadnavirus hominoidei",
                     "human coronavirus 229e" = "alphacoronavirus chicagoense",
                     "mammalian orthoreovirus" = "orthoreovirus mammalis",
                     "murine leukemia virus" = "gammaretrovirus murleu",
                     "venezuelan equine encephalitis virus" = "alphavirus venezuelan",
                     "western equine encephalitis virus" = "alphavirus western"
)

# then pull names
hvir <- humans$Virus

# filter bats by those names
bats %>% filter(Virus %in% hvir) -> zoo

# calculate richness
aggregate(Virus~Host,data=zoo, function(x) length(unique(x))) -> agg

# rename
agg <- agg %>% rename(zVirus = Virus)

# binary hosting
agg %>% mutate(dum_zvirus = 1) -> agg

# save 
write.csv(agg,"/Users/brianabetke/Desktop/bathaus/virus data/zoonotic viral response_bats VIRION flat.csv", row.names = FALSE)
