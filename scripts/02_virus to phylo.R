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
library(easyPubMed)

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
                          "Laephotis capensis"="Neoromicia capensis",
                          "Lissonycteris angolensis"="Myonycteris angolensis",
                          "Macronycteris commersoni"="Hipposideros commersoni",
                          "Macronycteris gigas"="Hipposideros gigas",
                          "Macronycteris vittata"="Hipposideros vittatus",
                          "Miniopterus africanus"="Miniopterus inflatus",
                          "Miniopterus orianae"="Miniopterus oceanensis",
                          "Molossus ater"="Molossus rufus",
                          "Myotis myotis/blythii"="Myotis myotis",
                          "Myotis oxygnathus"="Myotis blythii",
                          "Myotis ricketti"="Myotis pilosus",
                          "Neoromicia brunneus"="Neoromicia brunnea",
                          "Nyctalus velutinus"="Nyctalus plancyi",
                          "Parahypsugo crassulus"="Pipistrellus crassulus",
                          "Perimyotis subflavus"="Pipistrellus subflavus",
                          "Plecotus gaisleri"="Plecotus teneriffae",
                          "Pteronotus alitonus"="Pteronotus parnellii",
                          "Pteronotus rubiginosus"="Pteronotus parnellii",
                          "Rhinolophus blythi"="Rhinolophus lepidus",
                          #"Rhinolophus cornutus"="",
                          "Rhinolophus hildebrandti"="Rhinolophus hildebrandtii",
                          "Rhinolophus lobatus"="Rhinolophus landeri",
                          "Rhinolophus monoceros"="Rhinolophus pusillus",
                          "Rhinolophus rhodesiae"="Rhinolophus simulator"))

## recheck missing
setdiff(species$species,taxa$species)

## remove missing
species=species[!species$species%in%setdiff(species$species,taxa$species),]

## do we have duplicate species names after taxonomic matching?
table(table(species$species)>1)

## aggregate
species=aggregate(virus~species,data=species,sum)

## clean taxonomy file
taxa=taxa[c("species","gen","fam","clade","MSW3_sciName_matched","tip")]

## merge with species
data=merge(taxa,species,by="species",all.x=T)

## if no viruses, 0 viruses
data$virus=ifelse(is.na(data$virus),0,data$virus)

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

## merge
data=merge(data,cdata,by='tip')

## clean
rm(cites,vcites,x,cdata,i,counts)

## write flat file