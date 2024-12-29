# 10_IUCN data pull and merge.R
# Data pull for IUCN data by iucn and merge

# Script for data cleaning
rm(list=ls()) 
graphics.off()

# packages
library(tidyverse)
library(rredlist)

# Update synurbic data
## load in previous trait data
data <- read_csv("/Volumes/BETKE 2021/bathaus/flat files/master data and PanTHERIA.csv")

# view table updates
table(data$Synurbic)

# complete
table(data$Complete) # should be all yes

# # remove NAs in iucn2020 to get only combine names, then do the same matching you did for COMBINE
# cdata <- merged %>% drop_na(iucn2020_binomial)
sum(is.na(data$iucn2020_binomial)) # 18 NAs 
#data <- zoonotic %>% drop_na(iucn2020_binomial) 

# only run full data pull if needed. If not, just read in pulled dataset to merge with full virus data
pull = "no"
if(pull == "yes"){
  # unchanged COMBINE names before merging with upham
  trait_data <- read_csv("~/Desktop/Bats and Viruses/bathaus/COMBINE Datasets/COMBINE_archives/trait_data_reported.csv") %>%
    filter(order == "Chiroptera")
  
  ## rredlist scrape
  bapi="insertyourownkey"
  
  ## run IUCN function
  rdata=lapply(trait_data$iucn2020_binomial,function(x){
    res=rl_search(name=x,key=bapi)$result; 
    print(x)
    return(res)
  }
  )
  
  ## rbind
  rdata=do.call(rbind.data.frame,rdata)
  rdata$iucn2020_binomial=rdata$scientific_name
  
  # duplicates
  rdata$iucn2020_binomial[duplicated(rdata$iucn2020_binomial)]
  # [1] "Dermanura watsoni"      "Harpiocephalus harpia"  "Lonchophylla thomasi"  
  # [4] "Lophostoma carrikeri"   "Molossus coibensis"     "Myotis formosus"       
  # [7] "Myotis simus"           "Natalus mexicanus"      "Natalus stramineus"    
  # [10] "Nyctophilus corbeni"    "Pipistrellus kuhlii"    "Pteropus chrysoproctus"
  # [13] "Pteropus pelewensis"    "Rhinolophus borneensis" "Triaenops rufus" 
  
  # unique only
  rdata <- unique(rdata)
  
  # write as CSV
  write.csv(rdata, "/Users/brianabetke/Desktop/Bats and Viruses/bathaus/IUCN data/IUCN flat file.csv")

}else{ #read in the saved r data instead
# read in data instead if only need to merge
rdata <- read_csv("/Volumes/BETKE 2021/bathaus/IUCN data/IUCN flat file.csv")
}

# names 
# Matching the mismatches
rdata$mismatches <- rdata$iucn2020_binomial %>% 
  recode("Anoura carishina" = "Anoura canishina",
         "Casinycteris campomaanensis" = "Casinycteris ophiodon",
         "Chiroderma vizzotoi" = "Chiroderma vizottoi",
         "Dermanura azteca" = "Dermanura aztecus",
         "Dermanura cinerea" = "Dermanura cinereus",
         "Dermanura glauca" = "Dermanura glaucus",
         "Dermanura gnoma" = "Dermanura gnomus", 
         "Dermanura rosenbergi" = "Dermanura rosenbergii",
         "Dermanura tolteca" = "Dermanura toltecus",
         "Pteropus leucopterus" = "Desmalopex leucopterus",
         "Diclidurus isabella" = "Diclidurus isabellus",
         "Paremballonura atrata" = "Emballonura atrata",
         "Paremballonura tiavato" = "Emballonura tiavato",
         "Rhyneptesicus nasutus" = "Eptesicus nasutus",
         "Hypsugo affinis" = "Falsistrellus affinis",
         "Macronycteris commersoni" = "Hipposideros commersoni",
         "Macronycteris gigas" = "Hipposideros gigas",
         "Hipposideros pendlebury" = "Hipposideros pendelburyi",
         "Macronycteris thomensis" = "Hipposideros thomensis",
         "Macronycteris vittatus" = "Hipposideros vittatus",
         "Lonchophylla cadenai" = "Hsunycteris cadenai",
         "Lonchophylla pattoni" = "Hsunycteris pattoni",
         "Pipistrellus anchietae" = "Hypsugo anchietae",
         "Lyroderma lyra" = "Megaderma lyra",
         "Gardnerycteris crenulatum" = "Mimon crenulatum", 
         "Gardnerycteris koepckeae"  = "Mimon koepckeae",
         #"Miniopterus schreibersii" = "Miniopterus fuliginosus",
         "Miniopterus orianae" = "Miniopterus oceanensis",
         "Ozimops loriae" = "Mormopterus loriae",       
         "Murina lorelieae" = "Murina loreliae",          
         "Lissonycteris angolensis" = "Myonycteris angolensis",
         "Mormopterus kalinowskii" = "Nyctinomops kalinowskii", 
         "Hypsugo alaschanicus" = "Pipistrellus alaschanicus",
         "Hypsugo anthonyi" = "Pipistrellus anthonyi",   
         "Hypsugo arabicus" = "Pipistrellus arabicus",
         "Hypsugo ariel" = "Pipistrellus ariel",       
         "Hypsugo cadornae" = "Pipistrellus cadornae",  
         "Hypsugo eisentrauti" = "Pipistrellus eisentrauti",
         "Hypsugo joffrei" = "Pipistrellus joffrei",     
         "Hypsugo kitcheneri" = "Pipistrellus kitcheneri",
         "Hypsugo lophurus" = "Pipistrellus lophurus",   
         "Hypsugo macrotis" = "Pipistrellus macrotis",
         "Hypsugo musciculus" = "Pipistrellus musciculus",  
         "Hypsugo pulveratus" = "Pipistrellus pulveratus",
         "Hypsugo savii" = "Pipistrellus savii",       
         "Perimyotis subflavus" = "Pipistrellus subflavus",
         "Hypsugo vordermanni" = "Pipistrellus vordermanni",
         "Baeodon alleni" = "Rhogeessa alleni",
         "Baeodon gracilis" = "Rhogeessa gracilis",      
         "Boneia bidens" = "Rousettus bidens",
         "Austronomus australis" = "Tadarida australis",      
         "Chaerephon jobimena" = "Tadarida jobimena",
         "Austronomus kuboriensis" = "Tadarida kuboriensis",          
         "Vampyriscus bidens" = "Vampyressa bidens",
         "Vampyriscus brocki" = "Vampyressa brocki",
         "Vampyriscus nymphaea" = "Vampyressa nymphaea"       
  ) 

# look at them 
c <- rdata %>% filter(iucn2020_binomial != mismatches) %>% select(iucn2020_binomial, mismatches)
rm(c)

# cnames (recode names in our dataset)
data$icnames <- recode(data$species,
                        "Dermanura incomitatus" = "Dermanura watsoni",
                        "Harpiocephalus mordax" = "Harpiocephalus harpia",
                        "Hsunycteris thomasi" = "Lonchophylla thomasi",
                        "Lophostoma aequatorialis" = "Lophostoma occidentalis",
                        "Lophostoma yasuni" = "Lophostoma carrikeri",
                        "Miniopterus fuliginosus" ="Miniopterus schreibersii",
                        "Molossus barnesi" = "Molossus coibensis",
                        "Myotis flavus" = "Myotis formosus",
                        "Myotis midastactus" = "Myotis simus",
                        "Natalus lanatus" = "Natalus mexicanus",
                        "Natalus saturatus" = "Natalus stramineus",
                        "Nyctophilus timoriensis" = "Nyctophilus corbeni",
                        "Paracoelops megalotis" = "Hipposideros Pomona",
                        "Pipistrellus deserti" = "Pipistrellus kuhlii",
                        "Pteropus argentatus" = "Pteropus chrysoproctus",
                        "Pteropus yapensis" = "Pteropus pelewensis",
                        "Triaenops menamena" = "Triaenops rufus",
                        "Rhinolophus chaseni" = "Rhinolophus borneensis",
                        "Triaenops rufus" = "Triaenops persicus",
                        "Myotis aelleni" = "Myotis chiloensis",
                        "Myotis hajastanicus" = "Myotis aurascens"
)

# create cnames in cdata
rdata$icnames <- rdata$mismatches

# merge (should be 1287 total obs)
icname_merge <- merge(data, rdata[c("icnames","category","population_trend")], by = "icnames", all.x = TRUE) # %>%
#   select(!iucn2020_binomial.y) %>%
#   rename(iucn2020_binomial = iucn2020_binomial.x)

## missing
table(is.na(icname_merge$category), icname_merge$Synurbic)

# write as csv for merge file
write.csv(icname_merge, "/Volumes/BETKE 2021/bathaus/flat files/IUCN data merge to zoonotic.csv")