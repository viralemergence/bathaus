## 02_Zoonotic Virus Pull
## pulling all zoonotic viruses from VIRION
## danbeck@ou.edu
## Updated 07/22/24 by Briana Betke (trim to iso and PCR)

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(vroom)

## load VIRION
#setwd("~/Desktop/Virion")
virion=vroom("/Volumes/BETKE 2021/bathaus/Verena data/Virion.v0.2.1.csv.gz")

## trim to host and virus names NCBI resolved
virion=virion[virion$HostNCBIResolved==T & virion$VirusNCBIResolved==T,]

## to bats and humans
virion=virion[virion$HostOrder=="chiroptera" | virion$Host=="homo sapiens",]

## trim to PCR and/or Isolation Detection Methods
virion=virion[virion$DetectionMethod == "PCR/Sequencing" | virion$DetectionMethod == "Isolation/Observation",]

## remove missing
virion=virion[!is.na(virion$Host),]
virion=virion[!is.na(virion$Virus),]

## remove missing virus family
virion=virion[!is.na(virion$VirusFamily),]

## unique records
virion$unique=with(virion,paste(Host,Virus))
virion=virion[!duplicated(virion$unique),]
length(unique(virion$unique))

## tabulate
table(virion$HostOrder)

## number of hosts and viruses
length(unique(virion$Host))
length(unique(virion$Virus))

## prune to Chiroptera
bats=virion %>% 
  filter(HostOrder == "chiroptera")

## prune to humans
humans=virion%>% 
  filter(Host == "homo sapiens")

## unique viruses
hvir=unique(humans$Virus)
rm(humans)

## which viruses infect humans
bats$zoonotic=ifelse(bats$Virus%in%hvir,1,0)

# zoonotic virus richness 
species <- bats %>% 
  group_by(Host) %>%
  summarize(zvirus = sum(zoonotic)) %>%
  ungroup()

# write to flat 
write.csv(species,"/Volumes/BETKE 2021/bathaus/virus data/zoonotic viral response_bats VIRION flat.csv", row.names = FALSE)

