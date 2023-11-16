## 01_virus pull
## pulling all bat viruses from VIRION
## danbeck@ou.edu

## Updated to VIRION latest release (v0.2.1) on 11/14/23 by Briana Betke

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(vroom)

## load VIRION
setwd("~/Desktop/Virion")
virion=vroom("/Volumes/BETKE 2021/bathaus/Verena data/Virion.v0.2.1.csv.gz")

## prune to Chiroptera
bats=virion %>% 
  filter(HostOrder == "chiroptera")

## trim to host and virus names NCBI resolved
bats=bats[bats$HostNCBIResolved==T & bats$VirusNCBIResolved==T,]

## virus richness 
bats$virus=1
species=aggregate(virus~Host,data=bats,sum)

## write to flat file
setwd("/Volumes/BETKE 2021/bathaus/virus data")
write.csv(species,"viral response_bats VIRION flat.csv")
