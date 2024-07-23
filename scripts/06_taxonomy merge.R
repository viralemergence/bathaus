## 06_Taxonomy Merge
## Merging combined Santini with taxonomy, citations, and virus data

## ready workspace
gc()
rm(list=ls())
graphics.off()

## libraries
library(tidyverse)
library(plyr)

## merge taxonomy/viruses with citation
setwd("/Volumes/BETKE 2021/bathaus/flat files")
data=read.csv("bathaus_virus to phylo backbone.csv")
cites=read.csv("bathaus citations.csv")
data=merge(data,cites,by="tip")
rm(cites)

## clean
data$X.x=NULL
data$X.y=NULL

## merge santini
setwd("/Volumes/BETKE 2021/bathaus/Santini et al 2018")
sdata=read.csv("Santini_dwellers and visitors.csv")

## check names
setdiff(sdata$Species_ph,data$tip)

## fix
sdata$tip=revalue(sdata$Species_ph,
                  c("Neoromicia_nanus"="Neoromicia_nana"))
sdata=sdata[c("tip","Synurbic")]

## merge
data=merge(data,sdata,by="tip",all.x=T)
rm(sdata)

## virus == 0 and cites > 0
data$filter=ifelse(data$virus==0 & data$cites==0,"drop","keep")
data$vfilter=ifelse(data$virus==0 & data$vcites==0,"drop","keep")

## save whole dataset
setwd("/Volumes/BETKE 2021/bathaus/flat files")
write.csv(data,"master data_1287 species.csv")

## trim based on broad filter
set=data[data$filter=="keep",]
write.csv(data,"filter data_600 species.csv")
