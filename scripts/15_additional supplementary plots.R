# 15_predictions
# Variable importance of anthropogenic roosting and human population density 
# measures
# briana.a.betke-1@ou.edu

# clean environment
rm(list=ls())
graphics.off()

# test code for pulling out all the relative importance info from every model
library(tidyverse)
library(plotrix)
library(patchwork)
#library(gbm)
library(rstatix)
#library(ggridges)
#library(ggrepel)

# setwd
setwd("/Users/brianabetke/Desktop/bathaus")
# data <- readRDS("flat files/log cleaned dataset 30 cutoff.rds")

###################### BRT Results
# Read in rds files for model outputs from drive
# richness models
vfams_brts <- readRDS("flat files/virus family richness with brts.rds")
no_vfams_brts <- readRDS("flat files/virus family richness without brts.rds")
# zoonotic proportion models
zfams_brts <- readRDS("flat files/zoonotic family richness with brts.rds")
no_zfams_brts <- readRDS("flat files/zoonotic family richness without brts.rds")
# virus reservoir
vbinary_brts <- readRDS("flat files/binary virus with brts.rds")
no_vbinary_brts <- readRDS("flat files/binary virus without brts.rds")
# zoonotic virus reservoir
zbinary_brts <- readRDS("flat files/zoonotic binary with brts.rds")
no_zbinary_brts <- readRDS("flat files/zoonotic binary without brts.rds")

## Human density metrics
# Var sums function
human <- function(data_name){
  
  # pull relative importance
  vinf <- lapply(data_name,function(x) x$rinf)
  
  # bind with rbind
  data_vinf <- do.call(rbind,vinf)
  
  # tidy output
  df_name <- data_vinf %>%
    group_by(var) %>%
    summarize(avg = mean(rel.inf),
              se = std.error(rel.inf)) %>%
    ungroup()
  
  # rank isn't too important here, just the human density vars
  # mean population density
  # %tile
  # 
  hd <- c("log_X27.2_HuPopDen_Mean_n.km2","X27.4_HuPopDen_Change",
          "log_X27.1_HuPopDen_Min_n.km2", "log_X27.3_HuPopDen_5p_n.km2",
          "Synurbic")
  
  df_name %>% filter(var %in% hd) -> df_hd
  
  # add column of clean variable names
  df_hd$names <- df_hd$var 
  
  df_hd$names <- df_hd$names %>% recode("log_X27.2_HuPopDen_Mean_n.km2" = "Log Mean Human Density",
                                      "X27.4_HuPopDen_Change" = "Human Density Change",
                                      "log_X27.1_HuPopDen_Min_n.km2" = "Log Minimum Human Density",
                                      "log_X27.3_HuPopDen_5p_n.km2" = "Human density 5th %tile",
                                      "Synurbic" = "Anthropogenic Roosting")
  
  df_hd %>% arrange(avg) -> df_hd
  
}

# pull anthro vars from all models
human(vbinary_brts) -> vb_hd
human(zbinary_brts) -> zb_hd
human(vfams_brts) -> vf_hd
human(zfams_brts) -> zf_hd

# add outcome names
vb_hd$outcome <- "Virus Host"
zb_hd$outcome <- "Zoonotic Host"
vf_hd$outcome <- "Family Richness"
zf_hd$outcome <- "Zoonotic Family Richness"

# bind
hd <- rbind(vf_hd, zf_hd, vb_hd, zb_hd)

hd$color <- ifelse(hd$var == "Synurbic","red",NA)

# plot
png("/Users/brianabetke/Desktop/bathaus/figs/figure S10.png", width=9,height=5.5,units="in",res=600)
ggplot(hd, aes(x = reorder(names, desc(names)), y = avg, fill = color)) +
  geom_bar(stat = "identity",
           show.legend = FALSE) +
  geom_errorbar(aes(ymin = avg-se, ymax = avg+se),width = 0.2) +
  facet_wrap(~outcome) +
  xlab(" ") +
  ylab("Average Relative Importance (%)") +
  theme_bw() +
  coord_flip()
dev.off()

## Comparison of variable importance across outcomes for the 15 most important
# variables that predict anthropogenic roosting ability
anthrovar <- function(data_name){
  
  # pull relative importance
  vinf <- lapply(data_name,function(x) x$rinf)
  
  # bind with rbind
  data_vinf <- do.call(rbind,vinf)
  
  # tidy output
  df_name <- data_vinf %>%
    group_by(var) %>%
    summarize(avg = mean(rel.inf),
              se = std.error(rel.inf)) %>%
    ungroup()
  
  # filter to top 15 plus anthro roosting
  top <- c("log_X26.1_GR_Area_km2","habitat_breadth_n",
           "X28.1_Precip_Mean_mm","X30.1_AET_Mean_mm",
           "litter_size_n","X30.2_PET_Mean_mm",
           "dphy_plant","log_adult_body_length_mm",
           "det_fruit","log_adult_mass_g",
           "X26.2_GR_MaxLat_dd","category",
           "log_adult_forearm_length_mm", "log_X27.2_HuPopDen_Mean_n.km2",
           "X26.3_GR_MinLat_dd", "Synurbic")
  
  # might need to turn top into an df with a rank column 
  topdf <- tibble(var = top, rank = 1:16)
  
  merge(topdf, df_name, by = "var") -> df_rank
  
  # add column of clean variable names
  df_rank$names <- df_rank$var 
  
  df_rank$names <- df_rank$names %>% recode("log_X26.1_GR_Area_km2" = "Log Geographic Area",
                                      "habitat_breadth_n" = "Habitat Breadth",
                                      "X28.1_Precip_Mean_mm" = "Mean Precipitation",
                                      "X30.1_AET_Mean_mm" = "Mean AET",
                                      "litter_size_n" = "Litter Size",
                                      "X30.2_PET_Mean_mm" = "Mean PET",
                                      "dphy_plant" = "Diet Plants",
                                      "log_adult_body_length_mm" = "Log Adult Body Length",
                                      "det_fruit" = "Diet Fruit",
                                      "log_adult_mass_g" = "Log Adult Mass",
                                      "X26.2_GR_MaxLat_dd" = "Maximum Lattitude",
                                      "category" = "IUCN Status",
                                      "log_adult_forearm_length_mm" = "Log Adult Forearm Length",
                                      "log_X27.2_HuPopDen_Mean_n.km2" = "Mean Human Density",
                                      "X26.3_GR_MinLat_dd" = "Minimum Lattitude",
                                      "Synurbic" = "Anthropogenic Roosting")
  
  df_rank %>% arrange(rank) -> df_rank
}

# pull anthro vars from all models
anthrovar(vbinary_brts) -> vb_top
anthrovar(zbinary_brts) -> zb_top
anthrovar(vfams_brts) -> vf_top
anthrovar(zfams_brts) -> zf_top

# add outcome names
vb_top$outcome <- "Virus Host"
zb_top$outcome <- "Zoonotic Host"
vf_top$outcome <- "Family Richness"
zf_top$outcome <- "Zoonotic Family Richness"

# bind
anth <- rbind(vf_top, zf_top, vb_top, zb_top)

# factor levels for plot
clean <- c("Log Geographic Area","Habitat Breadth", "Mean Precipitation", "Mean AET",
           "Litter Size", "Mean PET", "Diet Plants","Log Adult Body Length","Diet Fruit",
           "Log Adult Mass","Maximum Lattitude","IUCN Status","Log Adult Forearm Length",
           "Mean Human Density","Minimum Lattitude","Anthropogenic Roosting")

anth$names <- factor(anth$names, levels = clean)

# # facet ggplot
# png("/Users/brianabetke/Desktop/bathaus/top anthro traits.png", width=9,height=5.5,units="in",res=600)
# ggplot(anth, aes(x = forcats::fct_rev(names), y = avg, fill = names)) +
#   geom_bar(stat = "identity",
#            show.legend = FALSE) +
#   geom_errorbar(aes(ymin = avg-se, ymax = avg+se),width = 0.2) +
#   facet_wrap(~outcome) +
#   xlab("Anthropogenic roosting trait profile") +
#   ylab("Average Relative Importance") +
#   #scale_x_discrete(labels = anth$names) +
#   theme_bw() +
#   coord_flip()
# dev.off()

# color by top traits and anthro just to be clear?
anth$color <- ifelse(anth$var == "Synurbic","red",NA)

png("/Users/brianabetke/Desktop/bathaus/figs/figure S12.png", width=9,height=5.5,units="in",res=600)
ggplot(anth, aes(x = forcats::fct_rev(names), y = avg, 
                 fill=color)) +
  geom_bar(stat = "identity",
           show.legend = FALSE) +
  geom_errorbar(aes(ymin = avg-se, ymax = avg+se),width = 0.2) +
  facet_wrap(~outcome) +
  xlab("Anthropogenic roosting trait profile") +
  ylab("Average Relative Importance") +
  #scale_fill_manual(name = var, ) +
  #scale_x_discrete(labels = anth$names) +
  theme_bw() +
  coord_flip()
dev.off()
 
# # Var sums function
# varsums <- function(data_name, anthro){
# 
#   # pull relative importance
#   vinf <- lapply(data_name,function(x) x$rinf)
#   
#   # bind with rbind
#   data_vinf <- do.call(rbind,vinf)
#   
#   # tidy output
#   df_name <- data_vinf %>%
#     group_by(var) %>%
#     summarize(avg = mean(rel.inf)) %>%
#     ungroup() %>%
#     arrange(desc(avg))
#   
# # if anthro is equal to yes, do the types like this
#   
#   if(anthro == "yes") {
#     
#     df_name <- df_name %>% 
#       mutate(type = ifelse(startsWith(var, "X"),"Geographic",
#                            ifelse(startsWith(var, "det"),"forage", "other")))
#     
#     # Create type variable
#     df_name <- mutate(df_name, type = ifelse(startsWith(var, "X"), "Geographic", 
#                                              ifelse(startsWith(var, "det"), "Forage", 
#                                                     ifelse(var %in% c("Palearctic", "Neotropical", "Afrotropical","Indomalayan","Nearctic",
#                                                                       "Oceanian","Australasian","glaciation", "habitat_breadth_n","altitude_breadth_m","disected_by_mountains",
#                                                                       "log_X26.1_GR_Area_km2", "log_X27.2_HuPopDen_Mean_n.km2","log_lower_elevation_m","upper_elevation_m",
#                                                                       "log_X27.1_HuPopDen_Min_n.km2","log_X27.3_HuPopDen_5p_n.km2", "island_dwelling"),"Geographic",
#                                                            ifelse(var %in% c("dphy_invertebrate","dphy_plant","dphy_vertebrate","trophic_level","foraging_stratum"),"Forage",
#                                                                   ifelse(startsWith(var, "fam"),"Phylogeny",
#                                                                          ifelse(var %in% c("category","population_trend"),"Conservation", 
#                                                                                 ifelse(var %in% c("log_cites","log_vcites"),"Citation Count", "Life History"))))))))
#     
#     df_name$type[df_name$var =="Synurbic"] <- "Anthropogenic Roosting"
#     
#   } else {
#     
#     df_name <- df_name %>% 
#       mutate(type = ifelse(startsWith(var, "X"),"Geographic",
#                            ifelse(startsWith(var, "det"),"forage", "other")))
#     
#     # Create type variable
#     df_name <- mutate(df_name, type = ifelse(startsWith(var, "X"), "Geographic", 
#                                              ifelse(startsWith(var, "det"), "Forage", 
#                                                     ifelse(var %in% c("Palearctic", "Neotropical", "Afrotropical","Indomalayan","Nearctic",
#                                                                       "Oceanian","Australasian","glaciation", "habitat_breadth_n","altitude_breadth_m","disected_by_mountains",
#                                                                       "log_X26.1_GR_Area_km2", "log_X27.2_HuPopDen_Mean_n.km2","log_lower_elevation_m","upper_elevation_m",
#                                                                       "log_X27.1_HuPopDen_Min_n.km2","log_X27.3_HuPopDen_5p_n.km2", "island_dwelling"),"Geographic",
#                                                            ifelse(var %in% c("dphy_invertebrate","dphy_plant","dphy_vertebrate","trophic_level","foraging_stratum"),"Forage",
#                                                                   ifelse(startsWith(var, "fam"),"Phylogeny",
#                                                                          ifelse(var %in% c("category","population_trend"),"Conservation", 
#                                                                                 ifelse(var %in% c("log_cites","log_vcites"),"Citation Count", "Life History"))))))))
#     
#   }
#   
#   # sum across categories 
#   df_name %>% 
#     group_by(type) %>%
#     summarise(sum = sum(avg),
#               se = std.error(avg)) %>%
#     ungroup() %>%
#     arrange(desc(sum)) -> df_name
#   
# }
# 
# # save data for plots
# vbin <- varsums(vbinary_brts, "yes")
# zbin <- varsums(zbinary_brts, "yes")
# vfam <- varsums(vfams_brts, "yes")
# zfam <- varsums(zfams_brts, "yes")
# 
# # create oucome variable for binding
# vbin$outcome <- "Virus Host"
# zbin$outcome <- "Zoonotic Host"
# vfam$outcome <- "Virus Family Richness"
# zfam$outcome <- "Zoonotic Family Richness"
# 
# # bind
# all <- rbind(vbin,zbin,vfam,zfam)
# 
# # plot
# ggplot(all, aes(x = reorder(type, sum), y = sum, fill = type)) +
#   geom_bar(stat = "identity",
#            show.legend = FALSE) +
#   geom_errorbar(aes(ymin = sum-se, ymax = sum+se),width = 0.2) +
#   facet_wrap(~outcome) +
#   xlab("Trait Category") +
#   ylab("Sum of Relative Importance") +
#   theme_bw() +
#   coord_flip()
# 
# 
# # save data for plots
# nvbin <- varsums(no_vbinary_brts, "no")
# nzbin <- varsums(no_zbinary_brts, "no")
# nvfam <- varsums(no_vfams_brts, "no")
# nzfam <- varsums(no_zfams_brts, "no")
# 
# 
# # create outcome variable for binding
# nvbin$outcome <- "Virus Host"
# nzbin$outcome <- "Zoonotic Host"
# nvfam$outcome <- "Virus Family Richness"
# nzfam$outcome <- "Zoonotic Family Richness"
# 
# # bind
# nall <- rbind(nvbin,nzbin,nvfam,nzfam)
# 
# # make the fig ya dig
# png("/Users/brianabetke/Desktop/bathaus/Sum of relative importance.png", width=8,height=5.5,units="in",res=600)
# ggplot(nall, aes(x = reorder(type, sum), y = sum, fill = type)) +
#   geom_bar(stat = "identity",
#            show.legend = FALSE) +
#   geom_errorbar(aes(ymin = sum-se, ymax = sum+se),width = 0.2) +
#   facet_wrap(~outcome) +
#   xlab("Trait Category") +
#   ylab("Sum of Relative Importance") +
#   theme_bw() +
#   coord_flip()
# dev.off()
# 
# # lets remove citation count?
# all %>% filter(!type == "Citation Count") %>%
# ggplot(aes(x = reorder(type, sum), y = sum, fill = type)) +
#   geom_bar(stat = "identity",
#            show.legend = FALSE) +
#   geom_errorbar(aes(ymin = sum-se, ymax = sum+se),width = 0.2) +
#   facet_wrap(~outcome) +
#   xlab("Trait Category") +
#   ylab("Sum of Relative Importance") +
#   theme_bw() +
#   coord_flip()
# 
# # lets remove citation count?
# all %>% filter(!type == "Citation Count") %>%
#   ggplot(aes(x = reorder(type, sum), y = sum, fill = type)) +
#   geom_bar(stat = "identity",
#            show.legend = FALSE) +
#   geom_errorbar(aes(ymin = sum-se, ymax = sum+se),width = 0.2) +
#   facet_wrap(~outcome) +
#   xlab("Trait Category") +
#   ylab("Sum of Relative Importance") +
#   theme_bw() +
#   coord_flip()
# 
# # lets remove citation count?
# nall %>% filter(!type == "Citation Count") %>%
#   ggplot(aes(x = reorder(type, sum), y = sum, fill = type)) +
#   geom_bar(stat = "identity",
#            show.legend = FALSE) +
#   geom_errorbar(aes(ymin = sum-se, ymax = sum+se),width = 0.2) +
#   facet_wrap(~outcome) +
#   xlab("Trait Category") +
#   ylab("Sum of Relative Importance") +
#   theme_bw() +
#   coord_flip()
# 
# # looking at differences for the categories
# # virus host 
# vbin %>% rename(wsum = sum) -> vb 
# nvbin %>% rename(wosum = sum) -> nvb
# 
# vbdifs <- merge(vb[c("type","wsum")], nvb[c("type","wosum","outcome")], by = "type", all.x = TRUE)
# vbdifs$difs <- vbdifs$wsum - vbdifs$wosum
# vbdifs
# 
# 
# # zoonotic virus host 
# zbin %>% rename(wsum = sum) -> zb 
# nzbin %>% rename(wosum = sum) -> nzb
# 
# zbdifs <- merge(zb[c("type","wsum")], nzb[c("type","wosum","outcome")], by = "type", all.x = TRUE)
# zbdifs$difs <- zbdifs$wsum - zbdifs$wosum
# zbdifs
# 
# 
# # richness 
# vfam %>% rename(wsum = sum) -> vf 
# nvfam %>% rename(wosum = sum) -> nvf
# 
# vfdifs <- merge(vf[c("type","wsum")], nvf[c("type","wosum","outcome")], by = "type", all.x = TRUE)
# vfdifs$difs <- vfdifs$wsum - vfdifs$wosum
# vfdifs
# 
# # zoonotic richness
# zfam %>% rename(wsum = sum) -> zf 
# nzfam %>% rename(wosum = sum) -> nzf
# 
# zfdifs <- merge(zf[c("type","wsum")], nzf[c("type","wosum","outcome")], by = "type", all.x = TRUE)
# zfdifs$difs <- zfdifs$wsum - zfdifs$wosum
# zfdifs