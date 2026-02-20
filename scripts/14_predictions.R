# 14_predictions
# Host model predictions and maps
# briana.a.betke-1@ou.edu
# updated 02/20/26

# clear environment
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)
library(PresenceAbsence)

###### Model predictions
# read in model prediction csv files
setwd("/Users/brianabetke/Desktop/bathaus")
zbin_apreds <- read.csv("flat files/zoonotic virus host predictions.csv")
vbin_apreds <- read.csv("flat files/virus host predictions.csv")

set.seed(12345) # seed

#### threshold predictions
# run threshold test and MSS. if yes, if not, skip to MSS method
### testing thresholds
### comparing MSS (3) to 85%, 90%, and 95% sensitivity (10)

## virus models with roosting
ts.p95 <- optimal.thresholds(data.frame(vbin_apreds[,c('species','dum_virus','with')]),
                           threshold = 10001,
                           opt.methods = c(3,4,7,8,10),
                           req.sens = 0.95,
                           na.rm = TRUE)

ts.p90 <- optimal.thresholds(data.frame(vbin_apreds[,c('species','dum_virus','with')]),
                             threshold = 10001,
                             opt.methods = c(3,4,7,8,10),
                             req.sens = 0.90,
                             na.rm = TRUE)

ts.p85 <- optimal.thresholds(data.frame(vbin_apreds[,c('species','dum_virus','with')]),
                             threshold = 10001,
                             opt.methods = c(3,4,7,8,10),
                             req.sens = 0.85,
                             na.rm = TRUE)
# function to sum 
cut.p95 <- function(x) {sum(vbin_apreds$with[vbin_apreds$dum_virus==0] > x)}
cut.p90 <- function(x) {sum(vbin_apreds$with[vbin_apreds$dum_virus==0] > x)}
cut.p85 <- function(x) {sum(vbin_apreds$with[vbin_apreds$dum_virus==0] > x)}

sapply(unlist(ts.p95[2]), cut.p95)
sapply(unlist(ts.p90[2]), cut.p90)
sapply(unlist(ts.p85[2]), cut.p85)
# sensitivity of 85% brings the number of species to the same as MSS

## Virus models without
nts.p95 <- optimal.thresholds(data.frame(vbin_apreds[,c('species','dum_virus','without')]),
                             threshold = 10001,
                             opt.methods = c(3,4,7,8,10),
                             req.sens = 0.95,
                             na.rm = TRUE)

nts.p90 <- optimal.thresholds(data.frame(vbin_apreds[,c('species','dum_virus','without')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.90,
                              na.rm = TRUE)

nts.p85 <- optimal.thresholds(data.frame(vbin_apreds[,c('species','dum_virus','without')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.85,
                              na.rm = TRUE)

ncut.p95 <- function(x) {sum(vbin_apreds$without[vbin_apreds$dum_virus==0] > x)}
ncut.p90 <- function(x) {sum(vbin_apreds$without[vbin_apreds$dum_virus==0] > x)}
ncut.p85 <- function(x) {sum(vbin_apreds$without[vbin_apreds$dum_virus==0] > x)}

sapply(unlist(nts.p95[2]), ncut.p95)
sapply(unlist(nts.p90[2]), ncut.p90)
sapply(unlist(nts.p85[2]), ncut.p85)
# so sensitivity of 85% brings the number of species to the same as MSS

### Zoonotic hosts
## virus models with roosting
zts.p95 <- optimal.thresholds(data.frame(zbin_apreds[,c('species','dum_zvirus','with')]),
                             threshold = 10001,
                             opt.methods = c(3,4,7,8,10),
                             req.sens = 0.95,
                             na.rm = TRUE)

zts.p90 <- optimal.thresholds(data.frame(zbin_apreds[,c('species','dum_zvirus','with')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.90,
                              na.rm = TRUE)

zts.p85 <- optimal.thresholds(data.frame(zbin_apreds[,c('species','dum_zvirus','with')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.85,
                              na.rm = TRUE)


zcut.p95 <- function(x) {sum(zbin_apreds$with[zbin_apreds$dum_zvirus==0] > x)}
zcut.p90 <- function(x) {sum(zbin_apreds$with[zbin_apreds$dum_zvirus==0] > x)}
zcut.p85 <- function(x) {sum(zbin_apreds$with[zbin_apreds$dum_zvirus==0] > x)}

sapply(unlist(zts.p95[2]), zcut.p95)
sapply(unlist(zts.p90[2]), zcut.p90)
sapply(unlist(zts.p85[2]), zcut.p85)
# sensitivity of 85% brings the number of species to under

## without
# 95%
nzts.p95 <- optimal.thresholds(data.frame(zbin_apreds[,c('species','dum_zvirus','without')]),
                              threshold = 10001,
                              opt.methods = c(3,4,7,8,10),
                              req.sens = 0.95,
                              na.rm = TRUE)

nzts.p90 <- optimal.thresholds(data.frame(zbin_apreds[,c('species','dum_zvirus','without')]),
                               threshold = 10001,
                               opt.methods = c(3,4,7,8,10),
                               req.sens = 0.90,
                               na.rm = TRUE)

nzts.p85 <- optimal.thresholds(data.frame(zbin_apreds[,c('species','dum_zvirus','without')]),
                               threshold = 10001,
                               opt.methods = c(3,4,7,8,10),
                               req.sens = 0.85,
                               na.rm = TRUE)

nzcut.p95 <- function(x) {sum(zbin_apreds$without[zbin_apreds$dum_zvirus==0] > x)}
nzcut.p90 <- function(x) {sum(zbin_apreds$without[zbin_apreds$dum_zvirus==0] > x)}
nzcut.p85 <- function(x) {sum(zbin_apreds$without[zbin_apreds$dum_zvirus==0] > x)}

sapply(unlist(nzts.p95[2]), nzcut.p95) 
sapply(unlist(nzts.p90[2]), nzcut.p90)
sapply(unlist(nzts.p85[2]), nzcut.p85)

# clean
rm(list = ls()[!ls() %in% c("vbin_apreds","zbin_apreds")])

### threshold - Moving forward with MSS
# virus host with 
t.vmod <- optimal.thresholds(data.frame(vbin_apreds[,c('species','dum_virus','with')]),
                             threshold = 10001,
                             opt.methods = 3,
                             req.sens = 0.95,
                             na.rm = TRUE)

# virus host without
t.nvmod <- optimal.thresholds(data.frame(vbin_apreds[,c('species','dum_virus','without')]),
                              threshold = 10001,
                              opt.methods = 3,
                              req.sens = 0.95,
                              na.rm = TRUE)
# zoonotic with 
t.zmod <- optimal.thresholds(data.frame(zbin_apreds[,c('species','dum_zvirus','with')]),
                             threshold = 10001,
                             opt.methods = 3,
                             req.sens = 0.95,
                             na.rm = TRUE)
# zoonotic without
t.nzmod <- optimal.thresholds(data.frame(zbin_apreds[,c('species','dum_zvirus','without')]),
                             threshold = 10001,
                             opt.methods = 3,
                             req.sens = 0.95,
                             na.rm = TRUE)

# binary results
vbin_apreds %>% mutate(bin_with = with > t.vmod$with,
                       bin_without = without > t.nvmod$without) -> pred

zbin_apreds %>% mutate(bin_with = with > t.zmod$with,
                       bin_without = without > t.nzmod$without) -> zpred

# novel hosts
# virus
table(pred$with[pred$dum_virus==0] > t.vmod$with)
table(pred$without[pred$dum_virus==0] > t.nvmod$without)

# zoonotic
table(zpred$with[zpred$dum_zvirus==0] > t.zmod$with)
table(zpred$without[zpred$dum_zvirus==0] > t.nzmod$without)

# Looking at overlap
pred %>% filter(dum_virus == 0 & with > t.vmod$with) %>% pull(species) -> vnovel
pred %>% filter(dum_virus == 0 & without > t.nvmod$without) %>% pull(species) -> n_vnovel
setdiff(vnovel, n_vnovel) # 6 additional species
setdiff(n_vnovel, vnovel) # 1 not shared 
intersect(vnovel, n_vnovel) -> vsame # 127 in common

zpred %>% filter(dum_zvirus == 0 & with >= t.zmod$with) %>% pull(species) -> znovel
zpred %>% filter(dum_zvirus == 0 & without >= t.nzmod$without) %>% pull(species) -> n_znovel
setdiff(znovel, n_znovel) # 2 species different
setdiff(n_znovel, znovel) # 4 species that aren't in with
intersect(znovel, n_znovel) -> zsame # 150 in common

# Bin knowns and novel. Spell out roosting status for figs
pred %>% 
  mutate(status_w = ifelse(dum_virus == 1, "known", ifelse(dum_virus == 0 & bin_with == 1,"novel", "cut")),
         status_wout = ifelse(dum_virus == 1, "known", ifelse(dum_virus == 0 & bin_without == 1,"novel", "cut")),
         roost = ifelse(Synurbic == 1, "anthropogenic roosting", "natural roosting")) -> pred

zpred %>% 
  mutate(status_w = ifelse(dum_zvirus == 1, "known", ifelse(dum_zvirus == 0 & bin_with == 1,"novel", "cut")),
         status_wout = ifelse(dum_zvirus == 1, "known", ifelse(dum_zvirus == 0 & bin_without == 1,"novel", "cut")),
         roost = ifelse(Synurbic == 1, "anthropogenic roosting", "natural roosting")) -> zpred

# how many are anthropogenic? 
# overall virus
filter(pred, bin_with == 1 & status_w == "novel" & roost == "anthropogenic roosting") %>% nrow() #83
filter(pred, bin_without == 1 & status_wout == "novel" & roost == "anthropogenic roosting") %>% nrow() #76

# zoonotic
filter(zpred, bin_with == 1 & status_w == "novel" & roost == "anthropogenic roosting") %>% nrow() # 95
filter(zpred, bin_without == 1 & status_wout == "novel" & roost == "anthropogenic roosting") %>% nrow() #95

# Family and biogeographical realm break downs
# need to read in trait data before family dummys
traits <- readRDS("flat files/synurbic and traits only.rds")

# novel zoonotic
pred %>% filter(status_w == "novel" | status_wout == "novel") -> vnov
zpred %>% filter(status_w == "novel" | status_wout == "novel") -> znov

# include citations to see if predicted species are also poorly sampled.
merge(traits[c("species", "fam", "biogeographical_realm", "category", "population_trend", "cites", "vcites")], vnov, by = "species") -> vnov
merge(traits[c("species", "fam", "biogeographical_realm", "category", "population_trend", "cites", "vcites")], znov, by = "species") -> znov

# family breakdown
vnov %>% filter(status_w == "novel") %>% count(fam) %>% arrange(-n)
vnov %>% filter(status_wout == "novel") %>% count(fam) %>% arrange(-n)

znov %>% filter(status_w == "novel") %>% count(fam) %>% arrange(-n)
znov %>% filter(status_wout == "novel") %>% count(fam) %>% arrange(-n)

# biogeographical realms
# breakdown
vnov %>% filter(status_w == "novel") %>% separate_rows(biogeographical_realm, sep = ", ") %>% count(biogeographical_realm) %>% arrange(-n)
znov %>% filter(status_w == "novel") %>% separate_rows(biogeographical_realm, sep = ", ") %>% count(biogeographical_realm) %>% arrange(-n)

# ciations and biogeographical ranges
vnov %>% filter(status_w == "novel") %>% separate_rows(biogeographical_realm, sep = ", ") %>% mutate(citations = ifelse(cites == 0, "none", ifelse(cites == 1, "one", "more"))) -> v_cites
znov %>% filter(status_w == "novel") %>% separate_rows(biogeographical_realm, sep = ", ") %>% mutate(citations = ifelse(cites == 0, "none", ifelse(cites == 1, "one", "more"))) -> z_cites

prop.table(table(v_cites$citations, v_cites$biogeographical_realm),2)
prop.table(table(z_cites$citations, z_cites$biogeographical_realm),2)

## format for supplemental tables
vnov %>% 
  mutate(commonality = ifelse(species %in% vsame, "common", ifelse(status_w == "cut","without roosting only","with roosting only"))) %>%
  select(species, fam, roost, biogeographical_realm, commonality, cites) -> t1 

t1$roost[is.na(t1$roost)] <- "missing"
t1$biogeographical_realm <- as.character(t1$biogeographical_realm)
t1$biogeographical_realm[is.na(t1$biogeographical_realm)] <- "no data"

#save
write.csv(t1, "flat files/Table S3 predictions.csv", row.names = FALSE)

# Zoonotic
znov %>% 
  mutate(commonality = ifelse(species %in% zsame, "common", ifelse(status_w == "cut","without roosting only","with roosting only"))) %>%
  select(species, fam, roost, biogeographical_realm, commonality, cites) -> t2 

t2$roost[is.na(t2$roost)] <- "missing"
t2$biogeographical_realm <- as.character(t2$biogeographical_realm)
t2$biogeographical_realm[is.na(t2$biogeographical_realm)] <- "no data"

write.csv(t2, "flat files/Table S4 predictions.csv", row.names = FALSE)

###### Make some maps!
# libraries
library(fasterize)
#library(rgdal)
library(raster)
library(sf)
library(geodata) # loads in terra package

# map of novel hosts for with and without anthro
# reading IUCN file
# ranges obtained from https://www.iucnredlist.org/resources/spatial-data-download
iucn <- st_read("MAMMALS_TERRESTRIAL_ONLY/MAMMALS_TERRESTRIAL_ONLY.shp") %>%
  filter(order_=="CHIROPTERA")

# look at the unique names
iucn_names <- data.frame(unique(iucn$sci_name))

# how many mismatches are there for the entire dataset?
(alt_miss=setdiff(pred$species,iucn_names$unique.iucn.sci_name.)) # 139 mismatches

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
setdiff(pred$species,iucn$sci_name)

### on to the actual maps
## make a blank raster from https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/
#setwd("/Volumes/BETKE 2021/bathaus/alt_2-5m_bil")
r <- disaggregate(raster("alt_2-5m_bil/alt.bil")*0,2)

# pull known hosts
pred %>% filter(status_w == "known" & roost == "anthropogenic roosting") %>% pull(species) -> kn_anth
pred %>% filter(status_w == "known" & roost == "natural roosting") %>% pull(species) -> kn_nat

# pull novel hosts 
pred %>% filter(status_w == "novel" & roost == "anthropogenic roosting") %>% pull(species) -> nov_anth_w
pred %>% filter(status_w == "novel" & roost == "natural roosting") %>% pull(species) -> nov_nat_w
pred %>% filter(status_wout == "novel" & roost == "anthropogenic roosting") %>% pull(species) -> nov_anth_wout

# all bats for scaling by bat diversity
pred %>% pull(species) -> all_bats

# filter by known hosts
iucn_kanth <- iucn[iucn$sci_name %in% kn_anth,] # known anthropogenic
iucn_knat <- iucn[iucn$sci_name %in% kn_nat,] # known natural

# Filter by novel hosts
iucn_anth_w <- iucn[iucn$sci_name %in% nov_anth_w,] 
iucn_nat_w <- iucn[iucn$sci_name %in% nov_nat_w,] 
iucn_anth_out <- iucn[iucn$sci_name %in% nov_anth_wout,]

# filter to bats in upham
iucn_all <- iucn[iucn$sci_name %in% all_bats,] # all bat species

# create raster layers
# knowns
map_kanth <- fasterize(iucn_kanth, r, fun="sum")
map_knat <- fasterize(iucn_knat, r, fun="sum")

# novels..
map_anth_w <- fasterize(iucn_anth_w, r, fun="sum")
map_nat_w <- fasterize(iucn_nat_w, r, fun="sum")
map_anth_wout <- fasterize(iucn_anth_out, r, fun="sum")

# base map
map_all <- fasterize(iucn_all, r, fun="sum")

# zero for continental area
fix <- function(x) {sum(x,r,na.rm=TRUE)+r}

map_kanth <- fix(map_kanth)
map_knat <- fix(map_knat)

map_anth_w <- fix(map_anth_w)
map_nat_w <-fix(map_nat_w)
map_anth_wout <- fix(map_anth_wout)

map_all <- fix(map_all)

# divide rasters for bat diversity
kanth_div <- map_kanth/map_all
knat_div <- map_knat/map_all

anth_div_w <- map_anth_w/map_all
nat_div_w <- map_nat_w/map_all
anth_div_wout <- map_anth_wout/map_all

# Fix continent again
kanth_div <- fix(kanth_div)
knat_div <- fix(knat_div)
anth_div_w <- fix(anth_div_w)
nat_div_w <- fix(nat_div_w)
anth_div_wout <- fix(anth_div_wout)

# clean 
rm(map_kanth, map_knat, map_anth_w, map_anth_wout, map_nat_w, r, map_all)

# raster stack
raster::stack(anth_div_w, anth_div_wout) -> comp
raster::stack(kanth_div, anth_div_w) -> anthro
raster::stack(knat_div, nat_div_w) -> natural

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

# map building function
range_maps <- function(map_name, pal, roost, known=NULL){
  
  if(roost == "anthro"){
  
    if(known == "known"){
      # make the plot with known and novel cols
      rasterVis::levelplot(map_name,  
                           margin = FALSE,
                           col.regions = pal,
                           layout = c(2,1),
                           #at = seq(0, 15, 1),
                           names.attr = c("Known Hosts", "Novel Hosts"),
                           alpha = 0.5, 
                           scales=list(draw=FALSE),
                           par.strip.text=list(cex=0.75),
                           xlab = NULL, ylab = list("Anthopogenic roosting", cex = 0.75),
                           maxpixels = 5e6)
      
      
    }else{
      # make the with without comparison plots
      rasterVis::levelplot(map_name,  
                           margin = FALSE,
                           col.regions = pal,
                           layout = c(2,1),
                           #at = seq(0, 15, 1),
                           names.attr = c("With Anthropogenic Roosting", "Without Anthropogenic Roosting"),
                           alpha = 0.5, 
                           scales=list(draw=FALSE),
                           par.strip.text=list(cex=0.75),
                           xlab = NULL, ylab = NULL,
                           maxpixels = 5e6)
      
    }
  
  }else{
    
    rasterVis::levelplot(map_name,  
                         margin = FALSE,
                         col.regions = pal,
                         layout = c(2,1),
                         #at = seq(0, 15, 1),
                         #names.attr = c("With Anthropogenic Roosting", "Without Anthropogenic Roosting"),
                         alpha = 0.5, 
                         scales=list(draw=FALSE),
                         par.strip.text=list(cex=0),
                         xlab = NULL, ylab = list("Natural roosting", cex = 0.75),
                         maxpixels = 5e6)
  }
}

# Knows and novel 
# raster stacks
range_maps(anthro, green, "anthro", "known") -> anth_only
range_maps(natural, purple, "nat") -> nat_only

# clean up plots 
anth_only$par.settings$layout.heights[
  c('bottom.padding',
    'top.padding',
    'key.sub.padding',
    'axis.xlab.padding',
    'key.axis.padding'
    #'between'
  ) ] <- -0.15

anth_only$aspect.fill <- TRUE

nat_only$par.settings$layout.heights[
  c( 'bottom.padding',
     'top.padding',
     'key.sub.padding',
     'axis.xlab.padding',
     'key.axis.padding'
  ) ] <- -0.15

nat_only$aspect.fill <- TRUE

# save 
png("figs/figure S11.png", width=7,height=3.5,units="in",res=600)
cowplot::plot_grid(anth_only, nat_only, ncol = 1)
dev.off()

## new comparison of with and without
range_maps(comp, green, "anthro", "comp") -> comp_only

# clean up plots 
comp_only$par.settings$layout.heights[
  c('bottom.padding',
    'top.padding',
    'key.sub.padding',
    'axis.xlab.padding',
    'key.axis.padding'
    #'between'
  ) ] <- -0.15

comp_only$par.settings$layout.widths[
  c('ylab.axis.padding')
] <- -1.25

comp_only$aspect.fill <- TRUE

# save
png("figs/figure 6.png", width=7,height=2,units="in",res=600)
comp_only
dev.off()
