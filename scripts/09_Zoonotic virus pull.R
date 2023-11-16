# Zoonotic Virus Pull
# pulling all zoonotic viruses from VIRION
# Updated 11/14/23

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(vroom)

## load VIRION
setwd("~/Desktop/Virion")
virion=vroom("/Volumes/BETKE 2021/bathaus/Verena data/Virion.v0.2.1.csv.gz")

## trim to host and virus names NCBI resolved
virion=virion[virion$HostNCBIResolved==T & virion$VirusNCBIResolved==T,]

## to bats and humans
virion=virion[virion$HostOrder=="chiroptera" | virion$Host=="homo sapiens",]

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
  summarize(zvirus = sum(zoonotic))

# then save a s csv. Send the updated script to 
setwd("/Volumes/BETKE 2021/bathaus/virus data")
write.csv(species,"zoonotic viral response_bats VIRION flat.csv", row.names = FALSE)

