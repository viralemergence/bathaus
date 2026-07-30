## 01_virus pull
## pulling all bat viruses from VIRION
## danbeck@ou.edu

# updated by Briana Betke on 07/30/26

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

# update host names to capitals to be consistent with PREDICT names and upham
virion$Host <- Hmisc::capitalize(virion$Host)

## Initial cleaning
# Pull chiroptera
bats <- virion %>% filter(HostOrder == "chiroptera"| HostOrder =="Chiroptera")
#bats <- virion %>% filter(HostOrder == "chiroptera")

# NCBI ratified hosts and viruses
bats <- bats %>% filter(HostNCBIResolved==T & VirusNCBIResolved==T)

# clean out NAs
bats <- bats %>% drop_na(Host, Virus)

# trim to iso and pcr only
bats <- bats %>% filter(DetectionMethod == "PCR/Sequencing" | DetectionMethod == "Isolation/Observation")

# look at how many unique assocations are ictv ratified
bats %>% mutate(ID = paste(Host, Virus)) -> bats 
bats %>% filter(!duplicated(ID)) -> bats

table(bats$ICTVRatified)
# FALSE  TRUE 
# 2730   703 

# remove observatons ending with sp. in host or virus name then look again at how many obs there are
bats %>% filter(!grepl("sp\\.", Host)) -> bats
bats %>% filter(!grepl("sp\\.", Virus)) -> bats

table(bats$ICTVRatified)
# FALSE  TRUE 
# 2072   677

# look at all non ictv ratified
bats %>% filter(!ICTVRatified == TRUE) -> nonrat
# how many unique viruses? 

length(unique(nonrat$Virus))
#[1] 1173

# save DF of nonratified names to be looked up
lookup <- data.frame(unique(nonrat$Virus))

# rename the column and organize alphabetically
lookup %>% arrange(unique.nonrat.Virus.) -> lookup

# write csv. 
csv = "no" # change to yes to save 
if(csv == "yes"){
  write.csv(lookup, "/Users/brianabetke/Desktop/bathaus/virus data/viruses to look up on ictv.csv")
}
rm(csv)

# read in the renames
renames <- read_csv("/Users/brianabetke/Desktop/bathaus/virus data/Virus reconciliation - renaming.csv")

# check for duplicates
table(duplicated(renames$`virion names`))
# FALSE 
# 1173 

# table
table(renames$status)

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

# trim to ictv ratified
bats %>% filter(ICTVRatified == T) -> bats

# check for duplicates one more time
table(duplicated(bats$ID))

# aggregate
aggregate(Virus~Host,data=bats, function(x) length(unique(x))) -> agg

# binary hosting
agg %>% mutate(dum_virus = 1) -> agg

# save 
write.csv(agg,"/Users/brianabetke/Desktop/bathaus/virus data/viral response_bats VIRION flat.csv", row.names = FALSE)
