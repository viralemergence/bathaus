# 12_data cleaning
# # Script for data cleaning
# babetke@utexas.edu

rm(list=ls()) 
graphics.off()

# Packages
library(tidyverse)
library(fastDummies) # dummy coding family

### Read in data
zoonotic <- read_csv("/Volumes/BETKE 2021/bathaus/flat files/IUCN data merge to zoonotic.csv")

## check for complete - finished data should be only Yes
table(zoonotic$Complete)

## checking that PanTHERIA data is updated from -999 to NA
table(is.na(zoonotic$X26.1_GR_Area_km2))

## names of bats that are left to revisit - should return nothing when revisits are done
zoonotic[zoonotic$Complete == "revisit",]$species

## remove the extinct bats - read in taxonomy
tax <- read_csv("/Volumes/BETKE 2021/bathaus/phylos/taxonomy_mamPhy_5911species.csv")

# How many bats are extinct?
names <- tax[tax$ord == "CHIROPTERA" & tax$extinct. == 1, ]$Species_Name
# [1] "Desmodus_draculae" "Pteropus_brunneus" "Pteropus_pilosus"  "Pteropus_subniger" "Pteropus_tokudae" 

# remove underscores
names <- sub("_", " ", names)

# remove the extinct bats
data <- filter(zoonotic, !species %in% names)
rm(zoonotic, tax, names)

# filter out the additional extinct bats in IUCN
names <- data[!is.na(data$category) & data$category == "EX", ]$species

# remove the extinct bats
data <- filter(data, !species %in% names)
rm(names)

# glimpse and see variables
glimpse(data)
colnames(data)

# clean out cols
data <- data %>% 
  select(-c("vfilter","filter","...1","pcnames", "ccnames", "icnames"))

#### Dummy Families
## make binary columns for genus
dums <- dummy_cols(data["fam"])

## unique, does not appear to have duplicates?
dums <- dums[!duplicated(dums$fam),]

# factor all variables
# the for loop turned everything to NA and gave 19 warnings
dums <- dums %>% 
  mutate(across(where(is.numeric), factor))

## merge - variables dont stay as factors when you merge
data <- merge(data,dums,by="fam",all.x=T)
rm(dums)

##### break up Biogeographical realms
temp <- data[c("species","biogeographical_realm")] 

temp_wide <- temp %>%
  separate_rows(biogeographical_realm, sep = ", ") %>%
  drop_na() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = biogeographical_realm, values_from = value, values_fill = 0) %>%
  mutate(across(where(is.numeric), factor))

data <- merge(data, temp_wide, by = "species", all = TRUE)
rm(temp, temp_wide)

#### Factor COMBINE and remove unnecessary columns
data <- data %>% # Synurbic and variables that are factors according to COMBINE
  mutate(across(c("Synurbic","hibernation_torpor","fossoriality","trophic_level",
                  "foraging_stratum","activity_cycle", "freshwater", 
                  "marine","terrestrial_non-volant", "terrestrial_volant","island_endemicity",
                  "disected_by_mountains", "glaciation", "biogeographical_realm", "category", "population_trend"), 
                factor)) %>% 
  select(-c("MSW3_sciName_matched"))

length(colnames(data)) # 109 columns

# save before trimming
saveRDS(data, "/Volumes/BETKE 2021/bathaus/flat files/synurbic and traits only.rds")

#### Clean out variables with low variance and coverage
# Variation
## mode function
mode.prop <- function(x) {
  ux <- unique(x[is.na(x)==FALSE])
  tab <- tabulate(match(na.omit(x), ux))
  max(tab)/length(x[is.na(x)==FALSE])
}

## assess variation across columns
vars <- data.frame(apply(data,2,function(x) mode.prop(x)),
                   apply(data,2,function(x) length(unique(x))))

## get names
vars$variables <- rownames(vars)
names(vars) <- c("var","uniq","column")
vars$var <- round(vars$var,2)

## if homogenous (100%)
vars$keep <- ifelse(vars$var<1,"keep","cut")
#vars$keep=ifelse(vars$column%in%c('hPCR','competence',"fam"),'keep',vars$keep)
vars <- vars[order(vars$keep),]

table(vars$keep)
# cut keep 
# 18   87

## trim
keeps <- vars[-which(vars$keep=="cut"),]$column

## drop if no variation
data <- data[keeps]
rm(keeps,vars)

## assess missing values
mval=data.frame(apply(data,2,function(x) length(x[!is.na(x)])/nrow(data)))

## get names
mval$variables=rownames(mval)
names(mval)=c("comp","column")
mval$comp=round(mval$comp,2)

# ggplot of coverage
png("/Volumes/BETKE 2021/bathaus/figs/figure S1.png", width=9.5,height=5.5,units="in",res=600)
mval %>% 
  filter(!column %in% c("clade", "gen", "tip", "species", "fam", "virus", "zvirus", "Complete")) %>%
  ggplot(aes(comp)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept=0.30,linetype=2,linewidth=0.5) +
  #geom_vline(xintercept=0.25,linetype=2,size=0.5) +
  theme_bw() +
  labs(y="frequency",
       x="trait coverage across bat species") +
  scale_x_continuous(labels = scales::percent)
dev.off()

mval$keep=ifelse(mval$comp>=0.30,"keep","cut")
table(mval$keep)
# cut keep 
# 18   73 

## order
mval=mval[order(mval$comp),]
keeps=mval[-which(mval$keep=="cut"),]$column

# cleaned up table of keeps
coverage_table <- mval %>% 
  filter(keep == "keep") %>%
  rename(Variable = column,
         Coverage = comp) %>%
  filter(!Variable %in% c("clade", "gen", "tip", "species", "fam", "virus", 
                          "zvirus", "Complete", "iucn2020_binomial","biogeographical_realm")) %>%
  select(-keep) %>%
  relocate(Coverage, .after = Variable)

rownames(coverage_table) <- NULL

# save as csv
write.csv(coverage_table, "/Volumes/BETKE 2021/bathaus/flat files/coverage_table.csv", row.names = FALSE)

## drop if not well represented
data=data[keeps]
rm(mval,keeps,coverage_table)

## Clean out remaining variables
data <- data %>% 
  select(-c("tip","gen","fam","clade","iucn2020_binomial","biogeographical_realm"))

colnames(data) # resulting in 66 variables total, 64 covariates

# reordering so species and response variables are in front
data %>% select(species, Synurbic, virus, zvirus, cites, everything()) -> data

# Save as RDS
saveRDS(data, "/Volumes/BETKE 2021/bathaus/flat files/cleaned dataset 30 cutoff.rds")

## version with transformed vars
# look into distribution of continuous variables
# pull numeric vars
num <- select(data, where(is.numeric))

# remove % diet variables
num <- select(num, !starts_with(c("det","dphy")))

# remove ordinal type variables
num <- num %>% select(!c(litters_per_year_n, litter_size_n, island_dwelling, habitat_breadth_n))

# histograms
Hmisc::hist.data.frame(num)

# log plus constant
log_data <- data
log_data$log_cites <- log1p(log_data$cites)
log_data$log_vcites <- log1p(log_data$vcites)
log_data$log_lower_elevation_m <- log1p(log_data$lower_elevation_m)
log_data$log_X26.1_GR_Area_km2 <- log1p(log_data$X26.1_GR_Area_km2)
log_data$log_X27.1_HuPopDen_Min_n.km2 <- log1p(log_data$X27.1_HuPopDen_Min_n.km2)
log_data$log_X27.2_HuPopDen_Mean_n.km2 <- log1p(log_data$X27.2_HuPopDen_Mean_n.km2)
log_data$log_X27.3_HuPopDen_5p_n.km2 <- log1p(log_data$X27.3_HuPopDen_5p_n.km2)

# log no constant
log_data$log_adult_body_length_mm <- log10(log_data$adult_body_length_mm)
log_data$log_adult_forearm_length_mm <- log10(log_data$adult_forearm_length_mm)
log_data$log_adult_mass_g <- log10(log_data$adult_mass_g)

log_data <- log_data %>%
  select(!c(cites, vcites, lower_elevation_m, X26.1_GR_Area_km2,
            X27.1_HuPopDen_Min_n.km2, X27.2_HuPopDen_Mean_n.km2, X27.3_HuPopDen_5p_n.km2,
            adult_body_length_mm, adult_forearm_length_mm, adult_mass_g))

# Save as RDS
saveRDS(log_data, "/Volumes/BETKE 2021/bathaus/flat files/log cleaned dataset 30 cutoff.rds")
