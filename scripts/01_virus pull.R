## 01_virus pull
## pulling all bat viruses from VIRION
## danbeck@ou.edu

## Updated to VIRION latest release (v0.2.1) on 11/14/23 by Briana Betke
## Updated code on 07/22/24 by Briana Betke 

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(vroom)

## load VIRION
# setwd("~/Desktop/Virion")
virion=vroom("/Volumes/BETKE 2021/bathaus/Verena data/Virion.v0.2.1.csv.gz")

## prune to Chiroptera
bats=virion %>% 
  filter(HostOrder == "chiroptera")

## trim to host and virus names NCBI resolved
bats %>% filter(HostNCBIResolved==T & VirusNCBIResolved==T) -> bats

## count the number of unique pairs per host species
vh_uni <- aggregate(Virus~Host,data=bats, function(x) length(unique(x)))
sum(vh_uni$Virus) #3775 total unique pairs.

## how many observations do not have viral families?
bats %>% drop_na(Host, Virus) %>% count(is.na(VirusFamily)) # 38

# Are there more than one detection per pair with missing virus family
bats %>% drop_na(Host, Virus) %>% mutate(ID = paste(Host, Virus)) %>% 
  group_by(ID) %>% tally(is.na(VirusFamily)) %>% filter(n > 0) %>% rename(missing = n) -> miss
# 34 unique pairs with missing virus families

# add non missing to be sure?
bats %>% drop_na(Host, Virus) %>% mutate(ID = paste(Host, Virus)) %>% 
  group_by(ID) %>% tally(!is.na(VirusFamily)) %>% rename(non = n) -> non

merge(miss, non, by = "ID") -> dat

# number of detections per pair (not unique detection method)
bats %>% drop_na(Host, Virus) %>% mutate(ID = paste(Host, Virus)) %>% 
  group_by(ID) %>% count() %>% rename(total = n) -> tot

# merge by ID
dat <- merge(dat, tot, by = "ID") 
# it looks like none of the pairs have a mixture of detentions with virus family and without

# remove
rm(vh_uni,miss,tot,non,dat)

## look at tally of detection methods
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

# A tibble: 1 × 4
#   sero   iso unknown   pcr
#   <dbl> <dbl>   <dbl> <dbl>
# 1   18    20     315   654

# clean
rm(det)

# clean to pairs that don't include missing virus families, antibody only, and unknown only dectections.
bats %>%
  drop_na(Host, Virus) %>%
  filter(DetectionMethod == "PCR/Sequencing" | DetectionMethod == "Isolation/Observation") %>%
  filter(!is.na(VirusFamily)) %>%
  group_by(Host) %>%
  summarize(length(unique(Virus))) %>%
  rename("virus" = "length(unique(Virus))") %>% 
  ungroup() -> species

# write to flat file
write.csv(species,"/Volumes/BETKE 2021/bathaus/virus data/viral response_bats VIRION flat.csv")
