# Host model predictions and maps
# babetke@utexas.edu

# clear environment
rm(list=ls()) 
graphics.off()

# data wranglin
library(tidyverse)

# Model predictions
# read in model prediction csv files
zres_apreds <- read.csv("/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus host predictions.csv")
vres_apreds <- read.csv("/Volumes/BETKE 2021/bathaus/flat files/virus host predictions.csv")

#### Get cut offs 
## Overall virus hosts
# how many predicted hosts with roosting model?
quantile(vres_apreds$with[vres_apreds$dum_virus==0], c(0.90, 0.95))
table(vres_apreds$with[vres_apreds$dum_virus==0] >= 0.85)
# FALSE  TRUE 
# 804    94 
table(vres_apreds$with[vres_apreds$dum_virus==0] >= 0.89)
# FALSE  TRUE 
# 863    35 

# how many predicted hosts with non-roost model
quantile(vres_apreds$without[vres_apreds$dum_virus==0], c(0.90, 0.95))
table(vres_apreds$without[vres_apreds$dum_virus==0] >= 0.85)
# FALSE  TRUE 
# 802    96 

table(vres_apreds$without[vres_apreds$dum_virus==0] >= 0.88)
# FALSE  TRUE 
# 851    47 

# how many species do they share?
vres_apreds %>% filter(dum_virus == 0 & with > 0.89) %>% pull(species) -> vnames_w

vres_apreds %>% filter(dum_virus == 0 & without > 0.88) %>% pull(species) -> vnames

intersect(vnames, vnames_w) # 35 species 
# [1] "Barbastella leucomelas"   "Enchisthenes hartii"      "Eptesicus gobiensis"     
# [4] "Eumops dabbenei"          "Eumops underwoodi"        "Hipposideros ater"       
# [7] "Idionycteris phyllotis"   "Kerivoula lanosa"         "Miniopterus aelleni"     
# [10] "Miniopterus brachytragos" "Miniopterus egeri"        "Miniopterus fraterculus" 
# [13] "Miniopterus majori"       "Miniopterus paululus"     "Miniopterus petersoni"   
# [16] "Miniopterus shortridgei"  "Murina tubinaris"         "Myotis bocagii"          
# [19] "Myotis montivagus"        "Myotis tricolor"          "Nyctalus montanus"       
# [22] "Nycticeinops schlieffeni" "Nyctinomops aurispinosus" "Peropteryx kappleri"     
# [25] "Promops centralis"        "Promops nasutus"          "Rhinolophus osgoodi"     
# [28] "Rhinolophus swinnyi"      "Rhinolophus yunanensis"   "Rousettus lanosus"       
# [31] "Scotomanes ornatus"       "Scotophilus robustus"     "Scotophilus viridis"     
# [34] "Tadarida aegyptiaca"      "Triaenops rufus" 

setdiff(vnames, vnames_w) # 12 species ided by without 
# [1] "Chaerephon bivittatus"  "Chiroderma salvini"     "Hypsugo anchietae"     
# [4] "Miniopterus griffithsi" "Miniopterus newtoni"    "Murina huttoni"        
# [7] "Phylloderma stenops"    "Pipistrellus rusticus"  "Pipistrellus tenuis"   
# [10] "Pygoderma bilabiatum"   "Rhinolophus marshalli"  "Tadarida fulminans"

# How many are anthropogenic roosting?
vres_apreds %>% filter(dum_virus == 0 & with >= 0.89) %>% filter(Synurbic == 1)
# 18 out of 35
vres_apreds %>% filter(dum_virus == 0 & without >= 0.88) %>% filter(Synurbic == 1)
# 23 of 47

## Zoonotic hosts
# with
quantile(zres_apreds$with[zres_apreds$dum_zvirus==0], c(0.90, 0.95))
table(zres_apreds$with[zres_apreds$dum_zvirus==0] > 0.71)
# FALSE  TRUE 
# 911   102 

table(zres_apreds$with[zres_apreds$dum_zvirus==0] > 0.77)
# FALSE  TRUE 
# 962    51 

# without
quantile(zres_apreds$without[zres_apreds$dum_zvirus==0], c(0.90, 0.95))
table(zres_apreds$without[zres_apreds$dum_zvirus==0] > 0.71)
# FALSE  TRUE 
# 910   103

table(zres_apreds$without[zres_apreds$dum_zvirus==0] > 0.77)
# FALSE  TRUE 
# 963    50 

# how many species do they share?
zres_apreds %>% filter(dum_zvirus == 0 & with > 0.77) %>% pull(species) -> znames_w

zres_apreds %>% filter(dum_zvirus == 0 & without > 0.77) %>% pull(species) -> znames

intersect(znames, znames_w) # 50 species so all that with predicted
# [1] "Barbastella leucomelas"    "Centurio senex"           
# [3] "Chiroderma salvini"        "Chiroderma villosum"      
# [5] "Cynomops paranus"          "Dermanura aztecus"        
# [7] "Enchisthenes hartii"       "Epomophorus crypturus"    
# [9] "Eumops dabbenei"           "Eumops underwoodi"        
# [11] "Glauconycteris variegata"  "Hipposideros ater"        
# [13] "Idionycteris phyllotis"    "Kerivoula lanosa"         
# [15] "Lonchophylla mordax"       "Macrophyllum macrophyllum"
# [17] "Murina huttoni"            "Myotis altarium"          
# [19] "Myotis bocagii"            "Myotis keaysi"            
# [21] "Myotis melanorhinus"       "Myotis montivagus"        
# [23] "Myotis oxyotus"            "Natalus lanatus"          
# [25] "Natalus mexicanus"         "Nycteris macrotis"        
# [27] "Nycticeinops schlieffeni"  "Nyctinomops aurispinosus" 
# [29] "Peropteryx kappleri"       "Peropteryx macrotis"      
# [31] "Phylloderma stenops"       "Phyllostomus elongatus"   
# [33] "Platyrrhinus vittatus"     "Promops centralis"        
# [35] "Promops nasutus"           "Pygoderma bilabiatum"     
# [37] "Rhinolophus fumigatus"     "Rhinolophus landeri"      
# [39] "Rhinolophus osgoodi"       "Rousettus lanosus"        
# [41] "Saccopteryx bilineata"     "Scotomanes ornatus"       
# [43] "Scotophilus dinganii"      "Scotophilus leucogaster"  
# [45] "Scotophilus viridis"       "Sturnira erythromos"      
# [47] "Tadarida aegyptiaca"       "Taphozous mauritianus"    
# [49] "Triaenops rufus"           "Uroderma magnirostrum"

setdiff(znames_w, znames) # one more from with than without
# "Myotis muricola" 

# How many are anthropogenic roosting?
zres_apreds %>% filter(dum_zvirus == 0 & with >= 0.77) %>% filter(Synurbic == 1)
# 30 out of 51 unknown
zres_apreds %>% filter(dum_zvirus == 0 & without >= 0.77) %>% filter(Synurbic == 1)
# 30 of 50 unknown

# Family and biogeographical realm break downs?
# need to read in trait data before dummys? 
traits <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/synurbic and traits only.rds")

# maybe match to these and get tables
vres_apreds %>% filter(dum_virus == 0 & with >= 0.89) -> vnovel

# include citations to see if predicted species are also poorly sampled.
merge(vnovel, traits[c("species", "fam", "biogeographical_realm", "cites", "vcites", "category")], by = "species") -> vnovel
table(vnovel$fam)
# EMBALLONURIDAE   HIPPOSIDERIDAE    MINIOPTERIDAE       MOLOSSIDAE   PHYLLOSTOMIDAE 
# 1                2                8                6                1 
# PTEROPODIDAE    RHINOLOPHIDAE VESPERTILIONIDAE 
# 1                3               13 

vnovel %>% separate_rows(biogeographical_realm, sep = ", ") -> vnovel
table(vnovel$biogeographical_realm)
# Afrotropical Australasian  Indomalayan     Nearctic  Neotropical   Palearctic 
# 15            2            8            3            7            5 

zres_apreds %>% filter(dum_zvirus == 0 & with >= 0.77) -> znovel

merge(znovel, traits[c("species", "fam", "biogeographical_realm")], by = "species") -> znovel
table(znovel$fam)
# EMBALLONURIDAE   HIPPOSIDERIDAE       MOLOSSIDAE        NATALIDAE       NYCTERIDAE 
# 4                2                7                2                1 
# PHYLLOSTOMIDAE     PTEROPODIDAE    RHINOLOPHIDAE VESPERTILIONIDAE 
# 13                2                3               17 

znovel %>% separate_rows(biogeographical_realm, sep = ", ") -> znovel
table(znovel$biogeographical_realm)
# Afrotropical  Indomalayan     Nearctic  Neotropical   Palearctic 
# 15            6            4           27            5 

#### maps
# you will need to get the binary status for known and unknown
vwith <- vres_apreds %>% 
  mutate(status = ifelse(dum_virus == 1, "known", ifelse(dum_virus == 0 & with > 0.89 ,"novel", "cut")),
         roost = ifelse(Synurbic == 1, "anthropogenic roosting", "natural roosting"))

vwith %>% filter(status != "cut") %>% drop_na(roost) %>% mutate(status = factor(status), roost = factor(roost))-> clean_vwith

table(vwith$status) # looks correct
# cut known novel 
# 863   381    35 

# read in bat shape files
bats=readRDS("/Volumes/BETKE 2021/bathaus/bat shp.rds")

## make species names match your dataset format
bats$tip=gsub("_"," ",bats$binomial)

## check missing
(miss=setdiff(clean_vwith$species,bats$tip))

# [1] "Dermanura cinereus"        "Dermanura glaucus"        
# [3] "Dermanura toltecus"        "Hipposideros commersoni"  
# [5] "Hipposideros gigas"        "Hipposideros vittatus"    
# [7] "Hsunycteris thomasi"       "Megaderma lyra"           
# [9] "Mimon crenulatum"          "Miniopterus fuliginosus"  
# [11] "Miniopterus mossambicus"   "Myonycteris angolensis"   
# [13] "Myotis flavus"             "Pipistrellus alaschanicus"
# [15] "Pipistrellus pulveratus"   "Pipistrellus savii"       
# [17] "Pipistrellus subflavus"    "Triaenops menamena"

# recode
bats$tip <- bats$tip %>% recode("Dermanura cinerea" = "Dermanura cinereus",
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
