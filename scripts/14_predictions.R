# Host model predictions and maps
# babetke@utexas.edu

# clear environment
rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)
library(PresenceAbsence)

# Model predictions
# read in model prediction csv files
zres_apreds <- read.csv("/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus host predictions.csv")
vres_apreds <- read.csv("/Volumes/BETKE 2021/bathaus/flat files/virus host predictions.csv")

set.seed(12345) # seed?

#### threshold predictions

### testing thresholds
# comparing MSS (3) to 85%, 90%, and 95% sensitivity (10)

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
# funtion to sum 
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

# virus host wihout
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

# how many are anthropogenic? 
pred %>% 
  mutate(status = ifelse(dum_virus == 1, "known", ifelse(dum_virus == 0 & bin_with == 1,"novel", "cut")),
         roost = ifelse(Synurbic == 1, "anthropogenic roosting", "natural roosting")) -> pred

table(pred$status, pred$roost, useNA = "ifany") # 65 anthropogenic

zpred %>% 
  mutate(status = ifelse(dum_zvirus == 1, "known", ifelse(dum_zvirus == 0 & bin_with == 1,"novel", "cut")),
         roost = ifelse(Synurbic == 1, "anthropogenic roosting", "natural roosting")) -> zpred

table(zpred$status, zpred$roost, useNA = "ifany") # 120 out of 189 thats almost 65%!!!

# # bar graph?
# zpred %>% filter(status == "novel") %>%
# ggplot(aes(x = roost, fill = roost)) +
#   geom_bar() +
#   scale_fill_manual(values = c("#8470ff","#9DD866","#A0B1BA")) +
#   theme_bw() +
#   theme(legend.position="none")

# Family and biogeographical realm break downs?
# need to read in trait data before dummys?
traits <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/synurbic and traits only.rds")

zpred %>% filter(status == "novel") -> novel

# props by roostin
prop.table(table(novel$roost, useNA = "ifany"))

# include citations to see if predicted species are also poorly sampled.
merge(novel, traits[c("species", "fam", "biogeographical_realm", "category")], by = "species") -> novel
table(novel$fam)
# EMBALLONURIDAE   HIPPOSIDERIDAE    MINIOPTERIDAE       MOLOSSIDAE   PHYLLOSTOMIDAE
# 1                2                8                6                1
# PTEROPODIDAE    RHINOLOPHIDAE VESPERTILIONIDAE
# 1                3               13

novel %>% separate_rows(biogeographical_realm, sep = ", ") %>% count(biogeographical_realm)

novel %>% filter(roost == "anthropogenic roosting") %>% separate_rows(biogeographical_realm, sep = ", ") -> anth_novel
ant_br <- data.frame(table(anth_novel$biogeographical_realm))
colnames(ant_br) <- c("realm","count")
ant_br$roost <- "anthropogenic roosting"

novel %>% filter(roost == "natural roosting") %>% separate_rows(biogeographical_realm, sep = ", ") -> nat_novel
nat_br <- data.frame(table(nat_novel$biogeographical_realm))
colnames(nat_br) <- c("realm","count")
nat_br$roost <- "natural roosting"

# bind
realms <- rbind(ant_br,nat_br)
ggplot(realms, aes(x = realm, fill = roost, weight = count, by = roost)) +
  geom_bar(position = "dodge") +
  scale_fill_manual(values = c("#8470ff","#9DD866","#A0B1BA")) +
  theme_bw() +
  theme(legend.position = "top")

#### maps
# you will need to get the binary status for known and unknown
# vwith <- vres_apreds %>% 
#   mutate(status = ifelse(dum_virus == 1, "known", ifelse(dum_virus == 0 & with > 0.89 ,"novel", "cut")),
#          roost = ifelse(Synurbic == 1, "anthropogenic roosting", "natural roosting"))

zpred %>% filter(status != "cut") %>% drop_na(roost) %>% mutate(status = factor(status), roost = factor(roost))-> clean

# table(vwith$status) # looks correct
# # cut known novel 
# # 863   381    35 

# read in bat shape files
bats=readRDS("/Volumes/BETKE 2021/bathaus/bat ranges/bat shp.rds")

## make species names match your dataset format
bats$tip=gsub("_"," ",bats$binomial)

## check missing
(miss=setdiff(clean$species,bats$tip))

# [1] "Dermanura cinereus"        "Dermanura glaucus"        
# [3] "Dermanura toltecus"        "Hipposideros commersoni"  
# [5] "Hipposideros gigas"        "Hipposideros vittatus"    
# [7] "Hsunycteris thomasi"       "Megaderma lyra"           
# [9] "Mimon crenulatum"          "Miniopterus fuliginosus"  
# [11] "Miniopterus mossambicus"   "Myonycteris angolensis"   
# [13] "Myotis flavus"             "Pipistrellus alaschanicus"
# [15] "Pipistrellus pulveratus"   "Pipistrellus savii"       
# [17] "Pipistrellus subflavus"    "Triaenops menamena"


# [1] "Dermanura glaucus"       "Dermanura toltecus"      "Harpiocephalus mordax"  
# [4] "Hipposideros commersoni" "Hipposideros gigas"      "Hipposideros vittatus"  
# [7] "Hsunycteris thomasi"     "Megaderma lyra"          "Mimon crenulatum"       
# [10] "Miniopterus fuliginosus" "Miniopterus mossambicus" "Myonycteris angolensis" 
# [13] "Myotis flavus"           "Myotis midastactus"      "Natalus lanatus"        
# [16] "Nyctophilus timoriensis" "Pipistrellus cadornae"   "Pipistrellus deserti"   
# [19] "Pipistrellus pulveratus" "Pipistrellus savii"      "Pipistrellus subflavus" 
# [22] "Pipistrellus tenuis"     "Triaenops menamena" 


# recode
bats$tip <- bats$tip %>% recode(#"Dermanura cinerea" = "Dermanura cinereus",
                    "Dermanura glauca" = "Dermanura glaucus",
                    "Dermanura tolteca" = "Dermanura toltecus",
                    "Macronycteris commersoni"= "Hipposideros commersoni",
                    "Macronycteris gigas" = "Hipposideros gigas",
                    "Macronycteris vittatus" = "Hipposideros vittatus",
                    "Lonchophylla thomasi" = "Hsunycteris thomasi",
                    "Lyroderma lyra" = "Megaderma lyra",
                    "Gardnerycteris crenulatum" = "Mimon crenulatum",
                    #"Miniopterus fuliginosus" = "Miniopterus schreibersii",
                    "Miniopterus orianae" = "Miniopterus oceanensis",
                    "Lissonycteris angolensis" = "Myonycteris angolensis",
                    "Hypsugo alaschanicus" = "Pipistrellus alaschanicus",
                    "Hypsugo pulveratus" = "Pipistrellus pulveratus",
                    "Hypsugo savii" = "Pipistrellus savii",
                    "Perimyotis subflavus" = "Pipistrellus subflavus")

## check missing
(miss=setdiff(clean_vwith$species,bats$tip))

## drop missing
clean_vwith=clean_vwith[!clean_vwith$species%in%miss,]

## trim
bats=bats[bats$tip%in%clean_vwith$species,]

## save id
bats$id=rownames(bats@data)

## simplify
library(rgeos)
tol=0.2

## loop through and simplify
lset=list()
for(i in 1:length(unique(bats$tip))){
  
  ## subset run
  set=bats[bats$tip==unique(bats$tip)[i],]
  
  ## save data
  sdata=set@data
  
  ## simplify
  shp=gSimplify(set,tol,topologyPreserve=TRUE)
  
  ## fortify
  shp=data.frame(fortify(shp,region="ID"))
  
  ## merge with sdata by id
  sdata=sdata[c("id","tip","binomial")]
  shp=merge(shp,sdata,by="id",all.x=T)
  
  ## save
  lset[[i]]=shp
  
  ## print
  print(paste(i,"in",length(unique(bats$tip))))
}

## convert to data
bset=do.call(rbind.data.frame,lset)

## clean
rm(bats)

## merge with data
clean_vwith %>% rename(tip = species) -> clean_vwith
bats=merge(bset,clean_vwith,by="tip",all.x=T)

## get world map
library(ggalt)
require(proj4)
library(ggthemes)
library(viridis)
library(mapproj)
wdata=map_data("world")
#wdata=wdata[-which(wdata$region=='Antarctica'),]

alpha <- ifelse(bats$status == "known", 0.20, 0.40)

## plot (this may take some time)
png("/Volumes/BETKE 2021/bathaus/figs/figure 6.png", width=6,height=4,units="in",res=600)
ggplot(wdata,aes(long,lat))+
  
  ## base layer
  geom_polygon(aes(group=group),
               fill="grey90",colour="grey90",size=0.2)+
  
  ## add shapefiles
  geom_polygon(data=bats,
               aes(group=paste(tip,group),
                   fill=roost), alpha = alpha) +
  facet_grid(roost ~ status) +
  #guides(fill="none") +
  theme_bw() +
  scale_fill_manual(labels = c("anthropogenic roosting","natural roosting"),
                    values = c("#9DD866","#8470ff")) +
  #scale_alpha_manual(values = c(0.20, 0.25)) +
  labs(fill = "Roosting") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position = "top",
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.box.spacing = unit(0, "cm")) +
  coord_map("mercator",xlim=c(-180,180)) 
dev.off()

# quick little bar plot for esa
zpred %>% select(bin_with, status, roost) %>% filter(status == "novel") -> hm
hm$model <- "with"
colnames(hm) <- c("bin", "status", "roost","model")

zpred %>% select(bin_without, status, roost) %>% filter(status == "novel" & bin_without == 1) -> hm2
hm2$model <- "without"
colnames(hm2) <- c("bin", "status", "roost","model")

# bind 
hm_tot <- rbind(hm,hm2)

# graph 
ggplot(hm_tot) +
  aes(x = model, fill = roost) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("#8470ff","#9DD866","#A0B1BA")) +
  theme_bw() +
  labs(y = "novel hosts", x = "models") #+
  #theme(legend.position = "top")

