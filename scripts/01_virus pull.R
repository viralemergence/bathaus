## 01_virus pull
## pulling all bat viruses from VIRION
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(vroom)

## load VIRION
setwd("~/Desktop/virion/Virion")
virion=vroom("Virion.csv.gz")

## prune to Chiroptera
bats=virion %>% 
  filter(HostOrder == "chiroptera")

## trim to host and virus names NCBI resolved
bats=bats[bats$HostNCBIResolved==T & bats$VirusNCBIResolved==T,]

## virus richness 
bats$virus=1
species=aggregate(virus~Host,data=bats,sum)

## write to flat file
setwd("~/Desktop/bathaus/virus data")
write.csv(species,"viral response_bats VIRION flat.csv")
