# 14_predictions
# Host model predictions and maps
# babetke@utexas.edu

# clear environment
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)
library(PresenceAbsence)

###### Model predictions
# read in model prediction csv files
zres_apreds <- read.csv("/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus host predictions.csv")
vres_apreds <- read.csv("/Volumes/BETKE 2021/bathaus/flat files/virus host predictions.csv")

set.seed(12345) # seed

#### threshold predictions
testrun = "no"
# run threshold test and MSS. if yes, if not, skip to MSS method
if(testrun == "yes"){
### testing thresholds
### comparing MSS (3) to 85%, 90%, and 95% sensitivity (10)

## virus models with roosting
ts.p95 <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','with')]),
                           threshold = 10001,
                           opt.methods = c(3,4,7,8,10),
                           req.sens = 0.95,
                           na.rm = TRUE)

ts.p90 <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','with')]),
                             threshold = 10001,
                             opt.methods = c(3,4,7,8,10),
                             req.sens = 0.90,
                             na.rm = TRUE)

ts.p85 <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','with')]),
                             threshold = 10001,
                             opt.methods = c(3,4,7,8,10),
                             req.sens = 0.85,
                             na.rm = TRUE)
# function to sum 
cut.p95 <- function(x) {sum(vres_apreds$with[vres_apreds$dum_virus==0] > x)}
cut.p90 <- function(x) {sum(vres_apreds$with[vres_apreds$dum_virus==0] > x)}
cut.p85 <- function(x) {sum(vres_apreds$with[vres_apreds$dum_virus==0] > x)}

sapply(unlist(ts.p95[2]), cut.p95)
sapply(unlist(ts.p90[2]), cut.p90)
sapply(unlist(ts.p85[2]), cut.p85)
# sensitivity of 85% brings the number of species to the same as MSS

## Virus models without
nts.p95 <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','without')]),
                             threshold = 10001,
                             opt.methods = c(3,4,7,8,10),
                             req.sens = 0.95,
                             na.rm = TRUE)

nts.p90 <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','without')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.90,
                              na.rm = TRUE)

nts.p85 <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','without')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.85,
                              na.rm = TRUE)

ncut.p95 <- function(x) {sum(vres_apreds$without[vres_apreds$dum_virus==0] > x)}
ncut.p90 <- function(x) {sum(vres_apreds$without[vres_apreds$dum_virus==0] > x)}
ncut.p85 <- function(x) {sum(vres_apreds$without[vres_apreds$dum_virus==0] > x)}

sapply(unlist(nts.p95[2]), ncut.p95)
sapply(unlist(nts.p90[2]), ncut.p90)
sapply(unlist(nts.p85[2]), ncut.p85)
# so sensitivity of 85% brings the number of species to the same as MSS

### Zoonotic hosts
## virus models with roosting
zts.p95 <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','with')]),
                             threshold = 10001,
                             opt.methods = c(3,4,7,8,10),
                             req.sens = 0.95,
                             na.rm = TRUE)

zts.p90 <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','with')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.90,
                              na.rm = TRUE)

zts.p85 <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','with')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.85,
                              na.rm = TRUE)


zcut.p95 <- function(x) {sum(zres_apreds$with[zres_apreds$dum_zvirus==0] > x)}
zcut.p90 <- function(x) {sum(zres_apreds$with[zres_apreds$dum_zvirus==0] > x)}
zcut.p85 <- function(x) {sum(zres_apreds$with[zres_apreds$dum_zvirus==0] > x)}

sapply(unlist(zts.p95[2]), zcut.p95)
sapply(unlist(zts.p90[2]), zcut.p90)
sapply(unlist(zts.p85[2]), zcut.p85)
# sensitivity of 85% brings the number of species to under

## without
# 95%
nzts.p95 <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','without')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.95,
                              na.rm = TRUE)

nzts.p90 <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','without')]),
                               threshold = 10001,
                               opt.methods = c(3,4,7,8,10),
                               req.sens = 0.90,
                               na.rm = TRUE)

nzts.p85 <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','without')]),
                               threshold = 10001,
                               opt.methods = c(3,4,7,8,10),
                               req.sens = 0.85,
                               na.rm = TRUE)

nzcut.p95 <- function(x) {sum(zres_apreds$without[zres_apreds$dum_zvirus==0] > x)}
nzcut.p90 <- function(x) {sum(zres_apreds$without[zres_apreds$dum_zvirus==0] > x)}
nzcut.p85 <- function(x) {sum(zres_apreds$without[zres_apreds$dum_zvirus==0] > x)}

sapply(unlist(nzts.p95[2]), nzcut.p95) 
sapply(unlist(nzts.p90[2]), nzcut.p90)
sapply(unlist(nzts.p85[2]), nzcut.p85)

# clean
rm(list = ls()[!ls() %in% c("vres_apreds","zres_apreds")])

### threshold - Moving forward with MSS
# virus host with 
t.vmod <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','with')]),
                             threshold = 10001,
                             opt.methods = 3,
                             req.sens = 0.95,
                             na.rm = TRUE)

# virus host without
t.nvmod <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','without')]),
                              threshold = 10001,
                              opt.methods = 3,
                              req.sens = 0.95,
                              na.rm = TRUE)
# zoonotic with 
t.zmod <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','with')]),
                             threshold = 10001,
                             opt.methods = 3,
                             req.sens = 0.95,
                             na.rm = TRUE)
# zoonotic without
t.nzmod <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','without')]),
                             threshold = 10001,
                             opt.methods = 3,
                             req.sens = 0.95,
                             na.rm = TRUE)
} else {
  ### threshold - Just MSS pull
  # virus host with 
  t.vmod <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','with')]),
                               threshold = 10001,
                               opt.methods = 3,
                               req.sens = 0.95,
                               na.rm = TRUE)
  
  # virus host without
  t.nvmod <- optimal.thresholds(data.frame(vres_apreds[,c('species','dum_virus','without')]),
                                threshold = 10001,
                                opt.methods = 3,
                                req.sens = 0.95,
                                na.rm = TRUE)
  # zoonotic with 
  t.zmod <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','with')]),
                               threshold = 10001,
                               opt.methods = 3,
                               req.sens = 0.95,
                               na.rm = TRUE)
  # zoonotic without
  t.nzmod <- optimal.thresholds(data.frame(zres_apreds[,c('species','dum_zvirus','without')]),
                                threshold = 10001,
                                opt.methods = 3,
                                req.sens = 0.95,
                                na.rm = TRUE)
}

# binary results
vres_apreds %>% mutate(bin_with = with > t.vmod$with,
                       bin_without = without > t.nvmod$without) -> pred

zres_apreds %>% mutate(bin_with = with > t.zmod$with,
                       bin_without = without > t.nzmod$without) -> zpred

# novel hosts
# virus
table(pred$with[pred$dum_virus==0] > t.vmod$with)
table(pred$without[pred$dum_virus==0] > t.nvmod$without)

# zoonotic
table(zpred$with[zpred$dum_zvirus==0] > t.zmod$with)
table(zpred$without[zpred$dum_zvirus==0] > t.nzmod$without)

# Looking at overlap
pred %>% filter(dum_virus == 0 & with >= t.vmod$with) %>% pull(species) -> vnovel
pred %>% filter(dum_virus == 0 & without >= t.nvmod$without) %>% pull(species) -> n_vnovel
setdiff(vnovel, n_vnovel) # 1 species different
intersect(vnovel, n_vnovel) # 110 in common

zpred %>% filter(dum_zvirus == 0 & with >= t.zmod$with) %>% pull(species) -> znovel
zpred %>% filter(dum_zvirus == 0 & without >= t.nzmod$without) %>% pull(species) -> n_znovel
setdiff(znovel, n_znovel) # 27 species different
intersect(znovel, n_znovel) # 162 in common

# Bin knowns and novel. Spell out roosting status for figs
pred %>% 
  mutate(status = ifelse(dum_virus == 1, "known", ifelse(dum_virus == 0 & bin_with == 1,"novel", "cut")),
         roost = ifelse(Synurbic == 1, "anthropogenic roosting", "natural roosting")) -> pred

zpred %>% 
  mutate(status = ifelse(dum_zvirus == 1, "known", ifelse(dum_zvirus == 0 & bin_with == 1,"novel", "cut")),
         roost = ifelse(Synurbic == 1, "anthropogenic roosting", "natural roosting")) -> zpred

# how many are anthropogenic? 
# overall virus
filter(pred, bin_with == 1 & status == "novel" & roost == "anthropogenic roosting") %>% pull(species) #65
filter(pred, bin_without == 1 & status == "novel" & roost == "anthropogenic roosting") %>% pull(species) #64

# zoonotic
filter(zpred, bin_with == 1 & status == "novel" & roost == "anthropogenic roosting") %>% pull(species) # 120
filter(zpred, bin_without == 1 & status == "novel" & roost == "anthropogenic roosting") %>% pull(species) #106

# Family and biogeographical realm break downs
# need to read in trait data before family dummys
traits <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/synurbic and traits only.rds")

# novel zoonotic
zpred %>% filter(status == "novel") -> novel

# props by roosting
prop.table(table(novel$roost, useNA = "ifany"))

# include citations to see if predicted species are also poorly sampled.
merge(traits[c("species", "fam", "biogeographical_realm", "category", "population_trend", "cites", "vcites")], novel, by = "species") -> novel

# family breakdown
table(novel$fam)

# Distribution of citations (all heavily right skewed)
hist(log(novel$cites+1))
hist(novel$vcites)

# biogeographical realms
novel %>% separate_rows(biogeographical_realm, sep = ", ") %>% count(biogeographical_realm)

# Conservation
table(novel$category)
table(novel$population_trend)

# conservation status v roosting ecology
table(novel$category, novel$roost, useNA = "ifany")
table(novel$population_trend, novel$roost, useNA = "ifany")

## look at traits of the 27 extra predicted by anth
names <- setdiff(znovel, n_znovel)

dif_traits <- traits %>% 
  filter(species %in% names)

# citation distributions
hist(dif_traits$cites)
hist(dif_traits$vcites)

# family
table(dif_traits$fam)

# roosting status
table(dif_traits$Synurbic) # fairly even distribution of nat (12) and anthro (14)

# geographic realm
dif_traits %>% separate_rows(biogeographical_realm, sep = ", ") %>% count(biogeographical_realm)

# conservation status
table(dif_traits$category)
table(dif_traits$population_trend)

###### Make some maps!!!
# libraries
library(fasterize)
library(rgdal)
library(raster)
library(sf)
library(geodata) # loads in terra package

# reading IUCN file
# ranges obtained from https://www.iucnredlist.org/resources/spatial-data-download
iucn <- st_read("/Volumes/BETKE 2021/bathaus/MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp") %>%
  filter(order_=="CHIROPTERA")

# look at the unique names
iucn_names <- data.frame(unique(iucn$sci_name))

# how many mismatches are there for the entire dataset?
(alt_miss=setdiff(zpred$species,iucn_names$unique.iucn.sci_name.)) # 139 mismatches

# reconcile for entire dataset
iucn$sci_name <- iucn$sci_name %>% # recode("old name" = "new name")
  recode("Thainycteris aureocollaris" = "Arielulus aureocollaris",
  "Mops aloysiisabaudiae" = "Chaerephon aloysiisabaudiae",
  "Mops ansorgei"= "Chaerephon ansorgei",
  "Mops atsinanana"  = "Chaerephon atsinanana",  
  "Mops bemmeleni" = "Chaerephon bemmeleni",
  "Mops bivittatus" = "Chaerephon bivittatus",
  "Mops bregullae" = "Chaerephon bregullae",
  "Mops chapini" = "Chaerephon chapini",
  "Mops gallagheri" = "Chaerephon gallagheri",
  "Mops jobensis" = "Chaerephon jobensis",
  "Mops johorensis" = "Chaerephon johorensis",
  "Mops major" = "Chaerephon major",
  "Mops nigeriae" = "Chaerephon nigeriae",
  "Mops plicatus" = "Chaerephon plicatus", 
  "Mops pumilus" = "Chaerephon pumilus",
  "Mops russatus" = "Chaerephon russatus",
  "Mops solomonis" = "Chaerephon solomonis",
  "Mops tomensis" = "Chaerephon tomensis",
  "Dermanura azteca" = "Dermanura aztecus",
  "Dermanura cinerea" = "Dermanura cinereus",
  "Dermanura glauca" = "Dermanura glaucus",
  "Dermanura gnoma" = "Dermanura gnomus",
  "Dermanura rosenbergi" = "Dermanura rosenbergii",
  "Dermanura tolteca" = "Dermanura toltecus",
  "Diaemus youngii" = "Diaemus youngi",
  "Diclidurus isabella" = "Diclidurus isabellus",
  "Paremballonura atrata" = "Emballonura atrata",
  "Paremballonura tiavato" = "Emballonura tiavato",
  "Epomophorus dobsonii" = "Epomops dobsonii",
  "Rhyneptesicus nasutus" = "Eptesicus nasutus",
  "Hypsugo affinis" = "Falsistrellus affinis",
  #"Glischropus aquilus"
  #"Harpiocephalus harpia" = "Harpiocephalus mordax",
  "Doryrhina camerunensis" = "Hipposideros camerunensis",
  "Macronycteris commersoni"= "Hipposideros commersoni",
  "Doryrhina cyclops" = "Hipposideros cyclops",
  "Macronycteris gigas" = "Hipposideros gigas",
  "Macronycteris thomensis" = "Hipposideros thomensis",
  "Macronycteris vittatus" = "Hipposideros vittatus",
  "Lonchophylla cadenai" = "Hsunycteris cadenai",
  "Lonchophylla pattoni" = "Hsunycteris pattoni",
  #"Lonchophylla thomasi" = "Hsunycteris thomasi",
  "Pipistrellus anchietae" = "Hypsugo anchietae",
  #"Lonchophylla inexpectata"
  "Lophostoma occidentalis" = "Lophostoma aequatorialis",
  # "Lophostoma yasuni"
  "Lyroderma lyra" = "Megaderma lyra",
  "Nesonycteris fardoulisi" = "Melonycteris fardoulisi",
  "Nesonycteris woodfordi" = "Melonycteris woodfordi",
  "Epomophorus intermedius" = "Micropteropus intermedius",
  "Epomophorus pusillus" = "Micropteropus pusillus",
  "Gardnerycteris crenulatum" = "Mimon crenulatum",
  "Gardnerycteris koepckeae"= "Mimon koepckeae",
  #"Miniopterus fuliginosus"
  #"Miniopterus mossambicus"
  "Miniopterus orianae" = "Miniopterus oceanensis",
  #"Molossus barnesi"
  "Ozimops beccarii" = "Mormopterus beccarii",
  "Setirostris eleryi" = "Mormopterus eleryi",
  "Ozimops halli" = "Mormopterus halli",
  "Ozimops kitcheneri" = "Mormopterus kitcheneri",
  "Ozimops loriae" = "Mormopterus loriae",
  "Ozimops lumsdenae" = "Mormopterus lumsdenae",
  "Micronomus norfolkensis" = "Mormopterus norfolkensis",
  "Ozimops planiceps" = "Mormopterus planiceps",
  "Murina cineracea" ="Murina feae",
  # Murina guilleni
  # Murina jaintiana
  "Murina lorelieae" = "Murina loreliae",
  #"Murina pluvialis"
  #"Murina tiensa"
  #"Myotis flavus",
  #"Myotis hajastanicus"
  #"Myotis handleyi"
  #"Myotis simus" = "Myotis midastactus",
  #"Myotis phanluongi"
  #"Natalus lanatus"
  "Notopteris macdonaldii" = "Notopteris macdonaldi",
  "Notopteris neocaledonicus" = "Notopteris neocaledonica",
  # Nycticeius aenobarbus
  "Nyctimene varius" = "Nyctimene minutus",
  "Nyctinomops kalinowskii" = "Mormopterus kalinowskii",
  #"Nyctophilus timoriensis"
  #"Paracoelops megalotis"
  "Paratriaenops furcula" = "Paratriaenops furculus",
  "Hypsugo alaschanicus" = "Pipistrellus alaschanicus",
  "Hypsugo anthonyi" = "Pipistrellus anthonyi",
  "Hypsugo arabicus" = "Pipistrellus arabicus",
  "Hypsugo ariel" = "Pipistrellus ariel",
  "Hypsugo cadornae" = "Pipistrellus cadornae",
  # "Pipistrellus deserti"
  "Hypsugo eisentrauti" = "Pipistrellus eisentrauti",
  "Hypsugo joffrei" = "Pipistrellus joffrei",
  "Hypsugo kitcheneri" = "Pipistrellus kitcheneri",
  "Hypsugo lophurus" = "Pipistrellus lophurus",
  "Hypsugo macrotis" = "Pipistrellus macrotis",
  "Hypsugo musciculus" = "Pipistrellus musciculus",
  "Hypsugo pulveratus" = "Pipistrellus pulveratus",
  "Hypsugo savii" = "Pipistrellus savii",
  "Perimyotis subflavus" = "Pipistrellus subflavus",
  #"Pipistrellus tenuis"
  "Hypsugo vordermanni"="Pipistrellus vordermanni",
  "Ptenochirus jagorii" = "Ptenochirus jagori",
  # "Pteropus argentatus"
  "Pteropus medius" = "Pteropus giganteus",
  "Pteropus vetula" = "Pteropus vetulus",
  #"Pteropus yapensis"
  #Rhinolophus chaseni
  #Rhinolophus francisi
  "Baeodon alleni" = "Rhogeessa alleni",
  "Baeodon gracilis" = "Rhogeessa gracilis",
  "Boneia bidens" = "Rousettus bidens",
  "Pilonycteris celebensis" = "Rousettus celebensis",
  "Austronomus australis" = "Tadarida australis",
  "Mops jobimena" = "Tadarida jobimena",
  "Austronomus kuboriensis" = "Tadarida kuboriensis",
  #"Triaenops parvus"
  #"Triaenops rufus"
  #"Uroderma bakeri"
  "Vampyriscus bidens" = "Vampyressa bidens",
  "Vampyriscus brocki" = "Vampyressa brocki",
  #"Vampyressa elisabethae"
  "Vampyriscus nymphaea" = "Vampyressa nymphaea",
  #Vampyressa sinchi
  )

# account for synonyms retained in the data
# pull names, change to synonym names and bind back
# list the names that need to be filtered and changed
cnames <- c("Harpiocephalus harpia","Lophostoma carrikeri","Lonchophylla thomasi","Lophostoma occidentalis",
            "Miniopterus schreibersii","Myotis formosus","Molossus coibensis","Murina harrisoni","Myotis aurascens",
            "Natalus mexicanus", "Myotis simus","Natalus stramineus","Nyctophilus corbeni","Hipposideros pomona",
            "Pipistrellus kuhlii","Pteropus chrysoproctus","Pteropus pelewensis","Rhinolophus borneensis","Triaenops persicus")

iucn_cnames <- iucn %>% 
  filter(sci_name %in% cnames)

iucn_cnames$sci_name <- iucn_cnames$sci_name %>% 
  recode("Harpiocephalus harpia" = "Harpiocephalus mordax",
  "Lophostoma carrikeri" = "Lophostoma yasuni",
  "Lonchophylla thomasi" = "Hsunycteris thomasi",
  "Lophostoma occidentalis" = "Lophostoma aequatorialis",
  "Miniopterus schreibersii" = "Miniopterus fuliginosus",
  "Myotis formosus" = "Myotis flavus",
  "Molossus coibensis" = "Molossus barnesi",
  "Murina harrisoni" = "Murina tiensa",
  "Myotis aurascens" = "Myotis hajastanicus",
  "Natalus mexicanus" = "Natalus lanatus", 
  "Myotis simus" = "Myotis midastactus",
  "Natalus stramineus" = "Natalus saturatus",
  "Nyctophilus corbeni" = "Nyctophilus timoriensis",
  "Hipposideros pomona" = "Paracoelops megalotis",
  "Pipistrellus kuhlii" = "Pipistrellus deserti",
  "Pteropus chrysoproctus" = "Pteropus argentatus",
  "Pteropus pelewensis"= "Pteropus yapensis",
  "Rhinolophus borneensis" = "Rhinolophus chaseni",
  "Triaenops persicus" = "Triaenops rufus"
)

# rbind to iucn
iucn <- rbind(iucn, iucn_cnames)

# check mismatches again - 33 
setdiff(zpred$species,iucn$sci_name)

### on to the actual maps
## make a blank raster from https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/
setwd("/Volumes/BETKE 2021/bathaus/alt_2-5m_bil")
r <- disaggregate(raster("alt.bil")*0,2)

# pull species names of novel hosts only
zpred %>% filter(status == "novel" & roost == "anthropogenic roosting") %>% pull(species) -> nov_anth
zpred %>% filter(status == "novel" & roost == "natural roosting") %>% pull(species) -> nov_nat
zpred %>% filter(status == "known" & roost == "anthropogenic roosting") %>% pull(species) -> kn_anth
zpred %>% filter(status == "known" & roost == "natural roosting") %>% pull(species) -> kn_nat
zpred %>% pull(species) -> all_bats

# filter by novel hosts
iucn_nanth <- iucn[iucn$sci_name %in% nov_anth,] # novel anthropogenic
iucn_nnat <- iucn[iucn$sci_name %in% nov_nat,] # novel natural
iucn_kanth <- iucn[iucn$sci_name %in% kn_anth,] # known anthropogenic
iucn_knat <- iucn[iucn$sci_name %in% kn_nat,] # known natural
iucn_all <- iucn[iucn$sci_name %in% all_bats,] # all bat species

# create raster layers
map_nanth <- fasterize(iucn_nanth, r, fun="sum")
map_nnat <- fasterize(iucn_nnat, r, fun="sum")
map_kanth <- fasterize(iucn_kanth, r, fun="sum")
map_knat <- fasterize(iucn_knat, r, fun="sum")
map_all <- fasterize(iucn_all, r, fun="sum")

# zero for continental area
fix <- function(x) {sum(x,r,na.rm=TRUE)+r}

map_nanth <- fix(map_nanth)
map_nnat <-fix(map_nnat)
map_kanth <- fix(map_kanth)
map_knat <- fix(map_knat)
map_all <- fix(map_all)

# divide rasters for bat diversity
nanth_div <- map_nanth/map_all
nnat_div <- map_nnat/map_all
kanth_div <- map_kanth/map_all
knat_div <- map_knat/map_all

# Fix continent again
nanth_div <- fix(nanth_div)
nnat_div <- fix(nnat_div)
kanth_div <- fix(kanth_div)
knat_div <- fix(knat_div)

# clean 
rm(map_nanth, map_nnat, map_kanth, map_knat, r, map_all)

# raster stack
raster::stack(kanth_div,nanth_div) -> anthro
raster::stack(knat_div, nnat_div) -> natural

# time to map!!
library(rasterVis)
library(RColorBrewer)

# color palettes
mycolors <- colorRampPalette(rev(brewer.pal(10,"Spectral")))(21)
mycolors[1] <- "#C0C0C0"

# anthro color palette 
green <- colorRampPalette(c("#9DD866", "#003300"))(15)
green <- append(green, "#C0C0C0", 0)
# view colors
scales::show_col(green)

# natural color palette
purple <- colorRampPalette(c("#907eff", "#0d0b19"))(15)
purple <- append(purple, "#C0C0C0", 0)
scales::show_col(purple)

# map function with custom color palettes for each plot
range_maps <- function(map_name, pal){
  
  rasterVis::levelplot(map_name,  
                       margin = FALSE,
                       col.regions = pal,
                       layout = c(2,1),
                       #at = seq(0, 15, 1),
                       alpha = 0.5, 
                       scales=list(draw=FALSE),
                       par.strip.text=list(cex=0),
                       xlab = NULL, ylab = NULL,
                       maxpixels = 5e6)
  
}

# raster stacks
range_maps(anthro, green) -> anth_only
range_maps(natural, purple) -> nat_only

# clean up plots 
anth_only$par.settings$layout.heights[
  c('bottom.padding',
    'top.padding',
    'key.sub.padding',
    'axis.xlab.padding',
    'key.axis.padding'
    #'between'
  ) ] <- -0.15
anth_only$par.settings$layout.widths[
  c('ylab.axis.padding')
] <- -1.25
anth_only$aspect.fill <- TRUE

nat_only$par.settings$layout.heights[
  c( 'bottom.padding',
     'top.padding',
     'key.sub.padding',
     'axis.xlab.padding',
     'key.axis.padding'
  ) ] <- -0.15
nat_only$par.settings$layout.widths[
  c('ylab.axis.padding')
] <- -1.25
nat_only$aspect.fill <- TRUE

# save
png("/Volumes/BETKE 2021/bathaus/figs/figure 5.png", width=7,height=3.5,units="in",res=600)
cowplot::plot_grid(anth_only, nat_only, ncol = 1)
dev.off()