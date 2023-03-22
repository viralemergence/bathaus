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
bats2=virion %>% 
  filter(HostOrder == "chiroptera")

## trim to host and virus names NCBI resolved
bats=bats[bats2$HostNCBIResolved==T & bats2$VirusNCBIResolved==T,]

## virus richness 
bats$virus=1
species2=aggregate(virus~Host,data=bats2,sum)

## write to flat file
setwd("~/Desktop/bathaus/virus data")
write.csv(species,"viral response_bats VIRION flat.csv")
