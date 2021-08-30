## 04_citation counts
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(Hmisc)
library(plyr)
library(easyPubMed)

## load in flat file
setwd("~/Desktop/bathaus/flat files")
data=read.csv("bathaus_virus to phylo backbone.csv")

## collect any citations per bat species
cites=c()
for(i in 1:length(data$tip)) {
  
  counts=as.numeric(as.character(get_pubmed_ids(gsub('_','-',data$tip[i]))$Count))
  cites[i]=counts
  print(paste(i,"/",nrow(data)))
}

## virus-related citations per bat species
vcites=c()
for(i in 1:length(data$tip)) {
  
  x=gsub('_','-',data$tip[i])
  x=paste("(",x,")",sep="")
  x=paste(x,"AND (virus OR viral)")
  counts=as.numeric(as.character(get_pubmed_ids(x)$Count))
  vcites[i]=counts
  print(paste(i,"/",nrow(data)))
}

## compile all citations
cdata=data.frame(tip=data$tip,
                 cites=cites,
                 vcites=vcites)

## clean
rm(cites,vcites,x,i,counts)

## write flat file
setwd("~/Desktop/bathaus/flat files")
write.csv(cdata,"bathaus citations.csv")