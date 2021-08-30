## 02_merge VIRION and Upham phylogeny
## danbeck@ou.edu

## clean environment & plots
rm(list=ls()) 
graphics.off()

## packages
library(tidyverse)
library(ape)
library(Hmisc)
library(plyr)

## load species-level virus data
setwd("~/Desktop/bathaus/virus data")
species=read.csv("viral response_bats VIRION flat.csv")

## load in Upham phylogeny
setwd("~/Desktop/bathaus/phylos")
tree=read.nexus('MamPhy_fullPosterior_BDvr_Completed_5911sp_topoCons_NDexp_MCC_v2_target.tre')

## load in taxonomy
taxa=read.csv('taxonomy_mamPhy_5911species.csv',header=T)
taxa=taxa[taxa$ord=="CHIROPTERA",]
taxa$tip=taxa$Species_Name

## trim phylo to bats
tree=keep.tip(tree,taxa$tiplabel)

## fix tip
tree$tip.label=sapply(strsplit(tree$tip.label,'_'),function(x) paste(x[1],x[2],sep=' '))
taxa$species=sapply(strsplit(taxa$tip,'_'),function(x) paste(x[1],x[2],sep=' '))

## correct VIRION names
species$species=capitalize(species$Host)
species$X=NULL
species$Host=NULL

## which names in VIRION aren't in Upham
setdiff(species$species,taxa$species)

## match manually
species$species=revalue(species$species,
                        c("Aeorestes cinereus"="Lasiurus cinereus",
                          "Aeorestes egregius"="Lasiurus egregius",
                          "Afronycteris nana"="Neoromicia nana",
                          "Antrozous dubiaquercus"="Bauerus dubiaquercus",
                          "Artibeus cinereus"="Dermanura cinereus",
                          "Artibeus glaucus"="Dermanura glaucus",
                          "Artibeus phaeotis"="Dermanura phaeotis",
                          "Artibeus toltecus"="Dermanura toltecus",
                          "Chaerephon leucogaster"="Chaerephon pumilus",
                          "Chaerephon pusillus"="Chaerephon pumilus",
                          "Dasypterus ega"="Lasiurus ega",
                          "Dasypterus intermedius"="Lasiurus intermedius",
                          "Dasypterus xanthinus"="Lasiurus xanthinus",
                          "Dobsonia andersoni"="Dobsonia anderseni",
                          "Dobsonia magna"="Dobsonia moluccensis",
                          "Doryrhina cyclops"="Hipposideros cyclops",
                          "Eptesicus regulus"="Vespadelus regulus",
                          "Eptesicus vulturnus"="Vespadelus vulturnus",
                          "Gardnerycteris crenulatum"="Mimon crenulatum",
                          "Hipposideros cf. ruber"="Hipposideros ruber",
                          "Hipposideros terasensis"="Hipposideros armiger",
                          "Hypsugo alaschanicus"="Pipistrellus alaschanicus",
                          "Hypsugo pulveratus"="Pipistrellus pulveratus",
                          "Hypsugo savii"="Pipistrellus savii",
                          "Laephotis capensis"="Eptesicus capensis",
                          "Lissonycteris angolensis"="Myonycteris angolensis",
                          "Macronycteris commersoni"=""
                          
                          
                          
                          
                          ))
