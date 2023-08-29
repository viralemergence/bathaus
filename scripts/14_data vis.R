# data vis - descriptive and brts?

#clean envrionment
rm(list=ls())
graphics.off()

# test code for pulling out all the relative importance info from every model
library(tidyverse)
library(plotrix)
library(patchwork)
library(gbm)
library(rstatix)
library(ggridges)
#library(ggrepel)

#### Descriptive stats
data <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/cleaned dataset 30 cutoff.rds")

# lab comp directory
data <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/cleaned dataset 30 cutoff.rds")

###################### BRT Results
# Read in rds files for model outputs from drive
vrichness_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/virus with brts.rds")
no_vrichness_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/virus without brts.rds")
zoo_prop_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/zoo_prop with brts.rds")
no_zoo_prop_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/zoo_prop without brts.rds")
vbinary_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/dum_virus with brts.rds")
no_vbinary_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/dum_virus without brts.rds")
zbinary_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/dum_zvirus with brts.rds")
no_zbinary_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/dum_zvirus without brts.rds")
cites_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/citation brts.rds")
vcites_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/virus citation brts.rds")

##### Summary stats of model performance
## Richness models - pseudo R2
# with
mean(sapply(vrichness_brts,function(x) x$testr2)) # 0.5575969
std.error(sapply(vrichness_brts,function(x) x$testr2)) # 0.0354378
# without
mean(sapply(no_vrichness_brts,function(x) x$testr2)) # 0.5544271
std.error(sapply(no_vrichness_brts,function(x) x$testr2)) # 0.0363453

## Zoonotic proportion - pseudo R2
# with
mean(sapply(zoo_prop_brts,function(x) x$testr2)) # 0.1333381
std.error(sapply(zoo_prop_brts,function(x) x$testr2)) # 0.00538303
# without
mean(sapply(no_zoo_prop_brts,function(x) x$testr2)) # 0.1331381
std.error(sapply(no_zoo_prop_brts,function(x) x$testr2)) # 0.005441455

## Virus reservoir - AUC
# With
mean(sapply(vbinary_brts,function(x) x$testAUC)) # 0.9030022
std.error(sapply(vbinary_brts,function(x) x$testAUC)) # 0.00176121
# without
mean(sapply(no_vbinary_brts,function(x) x$testAUC)) # 0.9029828
std.error(sapply(no_vbinary_brts,function(x) x$testAUC)) # 0.001763804

## Zoonotic virus reservoir - AUC
# With
mean(sapply(zbinary_brts,function(x) x$testAUC))
std.error(sapply(binary_brts,function(x) x$testAUC))
# without
mean(sapply(no_zbinary_brts,function(x) x$testAUC))
std.error(sapply(no_vbinary_brts,function(x) x$testAUC))

### ttests
## function for extracting data, perform unpaired t test, Cohen's d
# should probably take the with/without labels and model types
# adding data argument to allow for different responses
# Need to allow for removing negative pseudo R2 values
tfun=function(mod_with, mod_without, measure, fcol){
  
  ## format data
  n=length(sapply(mod_with,function(x) x$measure))
  adata=data.frame(y=c(sapply(mod_with,function(x) x[measure][[1]]),
                       sapply(mod_without,function(x) x[measure][[1]])),
                   response=c(rep('mod_with',n),rep('mod_without',n)),
                   seed=c(sapply(mod_with,function(x) x$seed),
                          sapply(mod_without,function(x) x$seed)))
  rm(n)
  
  # conditionally change negeative pseudo R2 values
  if(mod_with[[1]][["mod"]][["distribution"]][["name"]] != "bernoulli") {
    
    # change negatives to 0
    adata <- adata %>% mutate(y = ifelse(y <= 0, 0, y))
    
  }else{
    
    adata <- adata
    
  }
  
  ## factor
  adata$response=factor(adata$response,levels=c('mod_with','mod_without'))
  
  ## make jitter position
  adata$x=as.numeric(factor(adata$response))
  set.seed(1)
  adata$xj=jitter(adata$x,0.5)
  
  ## fix response
  adata$response2=recode(adata$response,
                         "mod_with"="with",
                         "mod_without"="without")
  
  ## t-test
  tsum=t.test(y~response,data=adata,
              alternative='two.sided',
              var.equal=F,paired=F)  
  
  ## effect size
  csum=cohens_d(y~response,data=adata,paired=F,var.equal=F)

  ## Add plot function?
  set.seed(3)
  plot <- ggplot(adata)+
    geom_boxplot(aes(x=x,y=y,group=x),width=0.25,outlier.alpha = 0,fill=fcol) +
    geom_point(aes(x=xj,y=y),size=1.5,alpha=0.5) +
    scale_x_continuous(breaks=c(1,2),
                       labels=levels(adata$response2),
                       limits=c(0.5,2.5)) +
    theme_bw() +
    theme(axis.text=element_text(size=10),
          axis.text.x=element_text(size=12),
          axis.title=element_text(size=12)) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    theme(axis.title.x=element_text(margin=margin(t=10,r=0,b=0,l=0))) +
    theme(axis.title.y=element_text(margin=margin(t=0,r=10,b=0,l=0))) +
    theme(legend.position = "none")
    # guides(colour="none")
  
  ## return
  return(list(adata=adata,tsum=tsum,csum=csum,plot=plot))
}

# virus richness models
vrichdata <- tfun(vrichness_brts, no_vrichness_brts, "testr2", "#E78AC3")

# view stats
vrichdata$tsum
vrichdata$csum

# boxplot with significance?
v_box <-vrichdata[["plot"]] + 
            labs(x = "Model Type", y = "Model Performance (Pseudo R2)", title = "Virus Richness") +
            geom_line(data = tibble(x=c(1, 2),y = c(0.97, 0.97)), aes(x=x,y=y), inherit.aes = FALSE) +
            geom_text(data = tibble(x=1.5,y = 0.995), 
                      aes(x=x,y=y, label = paste("T-test: p = ", round(vrichdata$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
            geom_text(data = tibble(x=1.5,y = 0.945), 
                      aes(x=x,y=y, label = paste("Cohen's d = ", round(vrichdata$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
            theme(plot.title = element_text(hjust = 0.5, size = 12)) 
v_box

rm(no_vrichness_brts)
rm(vrichness_brts)

## Zoonotic proportion
zpropdata <- tfun(zoo_prop_brts, no_zoo_prop_brts, "testr2", "#FC8D62")

# view stats
zpropdata$tsum
zpropdata$csum

# boxplot
z_box <- zpropdata[["plot"]] + labs(x = "Model Type", y = "Model Performance (Pseudo R2)") +
            labs(x = "Model Type", y = NULL, title = "Zoonotic Proportion Transformed") +
            geom_line(data = tibble(x=c(1, 2),y = c(0.29, 0.29)), aes(x=x,y=y), inherit.aes = FALSE) +
            geom_text(data = tibble(x=1.5,y = 0.3), 
                      aes(x=x,y=y, label = paste("T-test: p = ", round(zpropdata$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
            geom_text(data = tibble(x=1.5,y = 0.28), 
                      aes(x=x,y=y, label = paste("Cohen's d = ", round(zpropdata$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
            theme(plot.title = element_text(hjust = 0.5, size = 12))

z_box

rm(no_zoo_prop_brts, zoo_prop_brts)

v_box + z_box

## virus reservoir status
# need t stats for all metrics
vresAUC <- tfun(vbinary_brts, no_vbinary_brts, "testAUC", "#66C2A5")
vresSEN <- tfun(vbinary_brts, no_vbinary_brts, "sen")
vresSpec <- tfun(vbinary_brts, no_vbinary_brts, "spec")

# view stats
vresAUC$tsum
vresAUC$csum

# pvalue adjustment for values
ps=c(vresAUC$tsum$p.value,
     vresSEN$tsum$p.value,
     vresSpec$tsum$p.value)
round(p.adjust(ps,method="BH"),4)

# Plot AUC
# boxplots
vb_box  <- vresAUC[["plot"]] + labs(x = "Model Type", y = "Model Performance (Test AUC)", title = "Virus Host") +
              geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
              geom_text(data = tibble(x=1.5,y = 0.965), 
                        aes(x=x,y=y, label = paste("T-test: p = ", round(vresAUC$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
              geom_text(data = tibble(x=1.5,y = 0.955), 
                        aes(x=x,y=y, label = paste("Cohen's d = ", round(vresAUC$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
              theme(plot.title = element_text(hjust = 0.5, size = 12))
vb_box

rm(no_vbinary_brts, vbinary_brts)

# haven't changed the bar placement
vresSEN[["plot"]] + labs(x = "Model Type", y = "Model Performance (Sensitivity)") +
  geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.965), 
            aes(x=x,y=y, label = paste("T-test: p = ", round(vresSEN$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.955), 
            aes(x=x,y=y, label = paste("Cohen's d = ", round(vresSEN$csum$effsize, 4), sep = "")), inherit.aes = FALSE)

# haven't changed the bar placement
vresSpec[["plot"]] + labs(x = "Model Type", y = "Model Performance (Specificity)") +
  geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.965), 
            aes(x=x,y=y, label = paste("T-test: p = ", round(vresSpec$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.955), 
            aes(x=x,y=y, label = paste("Cohen's d = ", round(vresSpec$csum$effsize, 4), sep = "")), inherit.aes = FALSE)

## zoonotic reservoir status
# need t stats for all metrics
zresAUC <- tfun(zbinary_brts, no_zbinary_brts, "testAUC", "#8DA0CB")
zresSEN <- tfun(zbinary_brts, no_zbinary_brts, "sen")
zresSpec <- tfun(zbinary_brts, no_zbinary_brts, "spec")

# pvalue adjustment for values
ps=c(zresAUC$tsum$p.value,
     zresSEN$tsum$p.value,
     zresSpec$tsum$p.value)
round(p.adjust(ps,method="BH"),4)

# plot AUC
# boxplots
zresAUC[["plot"]] + labs(x = "Model Type", y = "Model Performance (Test AUC)") +
  geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.965), 
            aes(x=x,y=y, label = paste("T-test: p = ", round(vresAUC$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.955), 
            aes(x=x,y=y, label = paste("Cohen's d = ", round(vresAUC$csum$effsize, 4), sep = "")), inherit.aes = FALSE)

zresSEN[["plot"]] + labs(x = "Model Type", y = "Model Performance (Sensitivity)") + 
  geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.965), 
            aes(x=x,y=y, label = paste("T-test: p = ", round(vresAUC$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.955), 
            aes(x=x,y=y, label = paste("Cohen's d = ", round(vresAUC$csum$effsize, 4), sep = "")), inherit.aes = FALSE)

zresSpec[["plot"]] + labs(x = "Model Type", y = "Model Performance (Specificity)") +
  geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.965), 
            aes(x=x,y=y, label = paste("T-test: p = ", round(vresAUC$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
  geom_text(data = tibble(x=1.5,y = 0.955), 
            aes(x=x,y=y, label = paste("Cohen's d = ", round(vresAUC$csum$effsize, 4), sep = "")), inherit.aes = FALSE)

# patch them together and save
v_box + z_box / vb_box + zb_box

################### Variable Importance Plots and rankings
# Pull all the relative importance into a dataframe, get the mean, sd, and variation.
# Then create a plot similar to the one I made for the variants 
vinfPlot <- function(data_name, bar_color){
  
  # pull relative importance
  vinf <- lapply(data_name,function(x) x$rinf)
  
  # bind with rbind
  data_vinf <- do.call(rbind,vinf)
  
  # tidy output
  df_name <- data_vinf %>%
    group_by(var) %>%
    summarize(avg = mean(rel.inf),
              rse = std.error(rel.inf),
              rvar = var(rel.inf)) %>%
    ungroup() %>%
    arrange(desc(avg))
  
  #Clean up variable names
  df_name$var <- recode(df_name$var,
                         "cites" = "Citation Count",
                         "vcites" = "Virus Citation Count",
                         "X26.1_GR_Area_km2" = "Geographic Area",
                         "X30.1_AET_Mean_mm" = "Mean Monthly AET",
                         "X26.2_GR_MaxLat_dd" = "Maximum Latitude",
                         "X26.3_GR_MinLat_dd" = "Minimum Latitude",
                         "habitat_breadth_n" = "Habitat Breadth",
                         "litters_per_year_n" = "Litters Per Year",
                         "adult_body_length_mm" = "Adult Body Length",
                         "X28.2_Temp_Mean_01degC" = "Mean Monthly Temperature",
                         "X27.2_HuPopDen_Mean_n.km2" = "Mean Human Density",
                         "X28.1_Precip_Mean_mm" = "Mean Monthly Precipitation",
                         "litter_size_n" = "Litter Size",
                         "upper_elevation_m" = "Upper Elevation Limit",
                         "disected_by_mountains" = "Disected by Mountains",
                         "adult_forearm_length_mm" = "Adult Forearm Length",
                         "altitude_breadth_m" = "Altitude Breadth",
                         "X26.4_GR_MidRangeLat_dd" = "Median Latitudinal Range",
                         "foraging_stratum" = "Foraging stratum",
                         "adult_mass_g" = "Adult Mass",
                         "X30.2_PET_Mean_mm" = "Mean Monthly PET",
                         "det_vfish" = "Diet Fish",
                         "X26.5_GR_MaxLong_dd" = "Maximum Longitude",
                         "fam_RHINOLOPHIDAE" = "Rhinolophidae",
                         "det_diet_breadth_n" = "Diet Breadth",
                         "X26.6_GR_MinLong_dd" = "Minimum Longitude",
                         "det_vend" = "Diet Vend",
                         "X27.1_HuPopDen_Min_n.km2" = "Min Human Density",
                         "dphy_vertebrate" = "Diet Vertebrate",
                         "lower_elevation_m" = "Lower Elevation Limit",
                         "det_nect" = "Diet Nectar",
                         "X27.4_HuPopDen_Change" = "Human Density Change",
                         "X27.3_HuPopDen_5p_n.km2" = "Human Density 5th Percentile",
                         "X26.7_GR_MidRangeLong_dd" = "Median Longitudinal Range",
                         "det_fruit" = "Diet Fruit",
                         "fam_PHYLLOSTOMIDAE" = "Phyllostomidae",
                         "det_vect" = "Diet Vect",
                         "trophic_level" = "Trophic Level",
                         "dphy_invertebrate" = "Diet Invertebrate",
                         "island_dwelling" = "Island Dwelling",
                         "fam_MOLOSSIDAE" = "Molossidae",
                         "fam_MINIOPTERIDAE" = "Miniopteridae",
                         "fam_HIPPOSIDERIDAE" = "Hipposideridae",
                         "dphy_plant" = "Diet Plants",
                         "glaciation" = "Glaciation",
                         "fam_VESPERTILIONIDAE" = "Vespertilionidae",
                         "fam_EMBALLONURIDAE" = "Emballonuridae",
                         "fam_PTEROPODIDAE" = "Pteropodidae",
                         "activity_cycle" = "Activity Cycle",
                         "det_seed" = "Diet Seeds",
                         "fam_MORMOOPIDAE" = "Mormoopidae",
                         "fam_NATALIDAE" = "Natalidae",
                         "fam_NYCTERIDAE" = "Nycteridae",
                         "category" = "Conservation Status",
                         "population_trend" = "Population Trend",
                         "Synurbic" = "Anthropogenic Roost"
  )
  
  # you might want something that lets you also manipulate the color as well so you can 
  # make the points colored by dataset.
  # next thing would be to add the second dataset and see if I can have them on both
  
  # fig_name <- ggplot(df_name, aes(x = reorder(var, -avg), y = avg, 
  #                                       color = ifelse(var == "Synurbic", "Anthro", "Shade"))) + 
  #   #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
  #   geom_bar(stat = "identity") +
  #   geom_errorbar(aes(ymin = avg-rse, ymax = avg+rse))
  #   theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
  #         legend.position = "none") +
  #   scale_color_manual(breaks = c("Anthro","Shade"), values = c("red", box_color))
  
  fig_name <- ggplot(df_name, aes(x = reorder(var, avg), y = avg, 
                                  fill = ifelse(var == "Anthropogenic Roost", "Anthro", "Shade"))) + 
    #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = avg-rse, ymax = avg+rse)) +
    #geom_pointrange(aes(ymin = avg-rse, ymax = avg+rse)) +
    coord_flip() +
    theme_bw() +
    theme(# axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
          legend.position = "none") +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(size = 10)) + 
    scale_fill_manual(breaks = c("Anthro","Shade"), values = c("#F0027F", bar_color))
  
  # return a list with that dataset of rel.inf and figure
  return(list(df_name, fig_name))
}

# Get colors - try to avoid using the same colors as in Fig 1. Which doesn't leave a lot but maybe pastels 
library(RColorBrewer)
brewer.pal(n = 8, name = "Accent")

# run for all models with anthropogenic roosting
vrich_fig <- vinfPlot(vrichness_brts, "#E78AC3")
zprop_fig <- vinfPlot(zoo_prop_brts, "#FC8D62")
vbinary_fig <- vinfPlot(vbinary_brts, "#66C2A5")
zbinary_fig <- vinfPlot(zbinary_brts, "#8DA0CB")

# look at plots
vrich_gg <- vrich_fig[[2]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = "Relative Importance") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Virus Richness")

zprop_gg <- zprop_fig[[2]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = " ") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Zoonotic Proportion")

# Create the patchwork, dropping the y-axis labels from the plots, and setting
# the margins, this adds the common label
h_patch <- vrich_gg + zprop_gg & ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

# Use the tag label as a y-axis label
png("/Volumes/BETKE 2021/bathaus/figs/richnes_variableinf.png",width=7,height=6,units="in",res=600)
wrap_elements(h_patch) +
  labs(tag = "Relative Importance") +
  theme(
    plot.tag = element_text(size = 10),
    plot.tag.position = "bottom"
  )
dev.off()

# binary models
vbinary_gg <- vbinary_fig[[2]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = "Relative Importance") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Virus Host")

zbinary_gg <- zbinary_fig[[2]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = " ") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Zoonotic Host")

# Create the patchwork, dropping the y-axis labels from the plots, and setting
# the margins, this adds the common label
h_patch <- vbinary_gg + zbinary_gg & ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

# Use the tag label as a y-axis label
png("/Volumes/BETKE 2021/bathaus/figs/host_variableinf.png",width=7,height=6,units="in",res=600)
wrap_elements(h_patch) +
  labs(tag = "Relative Importance") +
  theme(
    plot.tag = element_text(size = 10),
    plot.tag.position = "bottom"
  )
dev.off()

## Rankings
# wouldn't you want to rank them from highest to lowest? 
virus <- vrich_fig[[1]]
virus$ranks <- 1:nrow(virus)
virus$type <- "virus richness"

zoop <- zprop_fig[[1]]
zoop$ranks <- 1:nrow(zoop)
zoop$type <- "zoonotic proportion"

vb <- vbinary_fig[[1]]
vb$ranks <- 1:nrow(vb)
vb$type <- "virus host"

zb <- zbinary_fig[[1]]
zb$ranks <- 1:nrow(zb)
zb$type <- "zoonotic host"

# merge into one dataset by variable
# try making lolipop plot for synurbic by filtering out anthropogenic roosting values and ggsegement - maybe match colors to var inf plot model colors?
ranks <- rbind(virus, zoop, vb, zb)
anthrank <- filter(ranks, var == "Anthropogenic Roost") %>% 
  mutate(group = ifelse(type == "virus richness" | type == "zoonotic proportion", "richness", "host"))

# geom segment
png("/Volumes/BETKE 2021/bathaus/figs/synurbic ranks.png",width=6.5,height=5,units="in",res=600)
ggplot(anthrank, aes(x=type, y=ranks, color = type)) +
  geom_segment(aes(x=type, xend=type, y=ranks, yend=0)) +
  geom_point(size=5) +
  coord_flip() +
  scale_y_reverse() +
  theme_bw() +
  labs(x = "Response", y = "Rank") +
  scale_color_manual(breaks = c("virus richness","zoonotic proportion","virus host","zoonotic host"), 
                     values = c("#E78AC3", "#FC8D62", "#66C2A5", "#8DA0CB")) +
  theme(panel.grid.major=element_blank())
dev.off()

# maybe add a smaller panel with citations? 


################# going to want code here that makes pdps for all models
# work on adapting the synurbat code to handle aggregated data - will aggregating result in a file that looks the same but just of averages?
# Need to look at hantavirus code 
# the thing is calculating the marginal effects across models do I actually show averaged y axis?


# so calculate pdep for every model
# Aggreagating pdps
pdp_agg=function(mod, feature){
  
  if(mod[["mod"]][["distribution"]][["name"]] == "gaussian"){
    
    pdep=plot.gbm(mod[["mod"]], 
                  i.var = feature,
                  return.grid = TRUE)
  }else{
  
  pdep=plot.gbm(mod[["mod"]], 
                i.var = feature, 
                type = "response",
                return.grid = TRUE)
  }
  
  pdep$seed=unique(mod[["seed"]])
  
  pdep$predictor = pdep[feature][,1]
  
  pdep$rank=1:nrow(pdep)
  
  pdep$yhat=pdep$y
  
  return(pdep)
  
}

## Partial dependence plots for continuous variables 
make_pdp_cont <- function(model,feature, pcolor) {
  
  agg = do.call(rbind, lapply(model,function(x) pdp_agg(x, feature)))
  
  ## get element-wise means
  x=with(agg,tapply(predictor,rank,mean))
  y=with(agg,tapply(yhat,rank,mean))
  
  ## save as mean
  pmean=data.frame(predictor=x,yhat=y)
  
  ## get yrange
  yrange=range(agg$yhat,pmean$yhat,na.rm=T)
  
  ## get histogram
  hi=hist(data[feature][,1],breaks=30,plot=F)
  hi=with(hi,data.frame(breaks[1:(length(breaks)-1)],counts))
  names(hi)=c("mids","counts")
  
  ## ggplot it
  ggplot(agg,aes(predictor,yhat,group=seed))+
    
    ## add histogram
    geom_segment(data=hi,inherit.aes=F,
                 aes(x=mids,xend=mids,
                     y=yrange[1],yend=plotrix::rescale(counts,yrange)),
                 size=5,colour="grey",alpha=0.50)+
    
    ## add lines
    geom_line(size=1,alpha=0.25,colour=pcolor)+
    
    ## add mean
    geom_line(data=pmean,size=2,inherit.aes=F,
              aes(predictor,yhat))+
    
    ## theme
    theme_bw()+
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=7))+
    theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0)))+
    theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0)))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
    labs(x=feature,y="marginal effect")+
    scale_y_continuous(labels=scales::number_format(accuracy=0.01))
  
}

## Function for factor pdp plots
make_pdp_fact <- function(model, feature, pcolor) {
  
  # aggregate
  agg = do.call(rbind, lapply(model,function(x) pdp_agg(x, feature)))
  
  ## get element-wise means
  y=with(agg,tapply(yhat,predictor,mean))
  
  ## save as mean
  #pmean=data.frame(predictor=x,yhat=y)
  pmean=data.frame(y)
  names(pmean)="yhat"
  pmean$predictor=rownames(pmean)
  rownames(pmean)=NULL
  
  ## make temp data
  temp=data
  temp$predictor=temp[feature][,1]
  
  ## do nothing
  agg=agg
  pmean=pmean
  temp=temp
  
  ## get yrange
  yrange=range(agg$yhat,pmean$yhat,na.rm=T)
  
  # pull counts for color
  df_cat <- as.data.frame(table(temp$predictor))
  
  # fix y axis point
  df_cat$ymin <- yrange[1]-0.01
  
  ## fix temp to yrange
  #temp$yhat=ifelse(temp$predictor==1,max(yrange),min(yrange))
  
  ## ggplot with rug
  set.seed(1)
  ggplot(agg,aes(predictor,yhat,group=seed)) +
    
    ## add individual BRTs
    geom_jitter(size=1,alpha=0.25,colour=pcolor,width=0.1) +
    
    ## add mean
    geom_point(data=pmean,size=2,inherit.aes=F,shape=15,
               aes(predictor,yhat)) +
    
    # # # ## add rug
    # geom_rug(data=df_cat,inherit.aes=F,
    #          aes(Var1, ymin),
    #          sides="b",position="jitter",
    #          colour="grey",alpha=0.25,
    #          na.rm=T)+
    
    geom_point(data=df_cat, inherit.aes = F, shape=23, 
               aes(x=Var1,y=ymin,fill=Freq)) +
    
    ## theme
    theme_bw() +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=7)) +
    theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0))) +
    theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0))) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    labs(x=feature,y="marginal effect") +
    scale_y_continuous(limits=c(yrange[1]-0.01,yrange[2]+0.01),
                       labels=scales::number_format(accuracy=0.01)) +
    scale_fill_continuous(high = "#525252", low = "#D9D9D9", guide="none")
  
}

# for right now, only grab synurbic pdps
vsyn <- make_pdp_fact(vrichness_brts, "Synurbic", "#E78AC3") 
zsyn <- make_pdp_fact(zoo_prop_brts, "Synurbic", "#FC8D62")
vbsyn <- make_pdp_fact(vbinary_brts, "Synurbic",  "#66C2A5")
zbsyn <- make_pdp_fact(zbinary_brts, "Synurbic",  "#8DA0CB")

png("/Volumes/BETKE 2021/bathaus/figs/synurbic pdps.png",width=6,height=4.5,units="in",res=600)
(vsyn + zsyn) / (vbsyn + zbsyn)
dev.off()

con <- make_pdp_fact(vrichness_brts, "category", "#E78AC3")
cit <- make_pdp_cont(vrichness_brts, "cites", "#E78AC3")
vit <- make_pdp_cont(vrichness_brts, "vcites", "#E78AC3")

## Plot partial dependence - Syunurbat code NEEDS TO BE ADAPTED
# No NA pdps
gr <- make_pdp_cont(noNA, "X26.1_GR_Area_km2", "Geographic Area (km2)", pcolor = FALSE)
hb <- make_pdp_cont(noNA,"habitat_breadth_n", "Habitat Breadth", pcolor = FALSE)
pm <- make_pdp_cont(noNA, "X28.1_Precip_Mean_mm", "Mean Monthly Precipitation (mm)", pcolor = FALSE)
at <- make_pdp_cont(noNA, "X30.1_AET_Mean_mm", "Mean Monthly AET", pcolor = FALSE)
ls <- make_pdp_cont(noNA, "litter_size_n", "Litter Size", pcolor = FALSE)
mp <- make_pdp_cont(noNA, "X30.2_PET_Mean_mm", "Mean Monthly PET", pcolor = FALSE)
am <- make_pdp_cont(noNA, "adult_mass_g","Adult Mass (g)", pcolor = FALSE)
dp <- make_pdp_cont(noNA, "dphy_plant","Diet Plants (%)", pcolor = FALSE)
bl <- make_pdp_cont(noNA, "adult_body_length_mm", "Adult Body Length", pcolor = FALSE)
fa <- make_pdp_cont(noNA, "adult_forearm_length_mm", "Adult Forearm Length", pcolor = FALSE)
mx <- make_pdp_cont(noNA, "X26.3_GR_MinLat_dd", "Maximum Latitude", pcolor = FALSE)
cs <- make_pdp_fact(noNA, "category", "Conservation Status", pcolor = FALSE)
fr <- make_pdp_cont(noNA, "det_fruit", "Diet Fruit (%)", pcolor = FALSE)
hp <- make_pdp_cont(noNA, "X27.2_HuPopDen_Mean_n.km2", "Mean Human Density", pcolor = FALSE)
ml <- make_pdp_cont(noNA, "X26.3_GR_MinLat_dd", "Minimum Latitude", pcolor = FALSE)

# Save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/Figure 3.png", width=7,height=7.5,units="in",res=300)
gr + hb + pm + at + ls + mp + am + dp + bl + fa + mx + cs + fr + hp + ml + plot_layout(nrow = 5, ncol = 3, byrow = TRUE)
dev.off()

png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/pdp 9 no NAs.png", width=10,height=8,units="in",res=300)
hb + gr + at + pm + ls + fr + fa + am + bl + plot_layout(nrow = 3, ncol = 3, byrow = TRUE)
dev.off()

#  Pseudo model pdps 
phb <- make_pdp_cont(pseudo,"habitat_breadth_n", "Habitat Breadth", pcolor = TRUE)
pgr <- make_pdp_cont(pseudo, "X26.1_GR_Area_km2", "Geographic Area (km2)", pcolor = TRUE)
ppm <- make_pdp_cont(pseudo, "X28.1_Precip_Mean_mm", "Mean Monthly Precipitation (mm)", pcolor = TRUE)
pat <- make_pdp_cont(pseudo, "X30.1_AET_Mean_mm", "Mean Monthly AET", pcolor = TRUE)
pls <- make_pdp_cont(pseudo, "litter_size_n", "Litter Size", pcolor = TRUE)
pmx <- make_pdp_cont(pseudo, "X26.3_GR_MinLat_dd", "Maximum Latitude", pcolor = TRUE)
pcs <- make_pdp_fact(pseudo, "category", "Conservation Status", pcolor = TRUE)
pcc <- make_pdp_cont(pseudo, "cites","Citation Count", pcolor = TRUE)
pam <- make_pdp_cont(pseudo, "adult_mass_g","Adult Mass (g)", pcolor = TRUE)
ppt <- make_pdp_cont(pseudo, "X30.2_PET_Mean_mm", "Mean Monthly PET", pcolor = TRUE)
pdp <- make_pdp_cont(pseudo, "dphy_plant","Diet Plants (%)", pcolor = TRUE)
pfr <- make_pdp_cont(pseudo, "det_fruit", "Diet Fruit (%)", pcolor = TRUE)
pml <- make_pdp_cont(pseudo, "X26.3_GR_MinLat_dd", "Minimum Latitude", pcolor = TRUE)
pfa <- make_pdp_cont(pseudo, "adult_forearm_length_mm", "Adult Forearm Length", pcolor = TRUE)

# Save
png("/Users/brianabetke/Desktop/Synurbic_Bats/synurbat/figures/Figure S2.png", width=7,height=7.5,units="in",res=300)
phb + pgr + pat + ppm + pat + pls + pmx + pcs + pcc + pam + ppt + pdp + pfr + pml + pfa + plot_layout(nrow = 5, ncol = 3, byrow = TRUE)
dev.off()

############## model predictions
# pulling predictions (may just want the values separately for now?)
# csv of predicted reservoirs and another for zoonotic?

## average predictions: Overall virus reservoir
virus_apreds=lapply(dum_virus_brts,function(x) x$predict)
virus_apreds=do.call(rbind,virus_apreds)

## aggregate
virus_apreds=data.frame(aggregate(pred~species,data=virus_apreds,mean),
                        # aggregate(cpred~treename,data=pcr_apreds,mean)['cpred'], ## holding wos constant
                        aggregate(dum_virus~species,data=virus_apreds,prod)["dum_virus"],
                        aggregate(dum_zvirus~species,data=virus_apreds,prod)["dum_zvirus"])

### type
# virus_apreds$type='PCR'

## average predictions: Zoonotic
zvirus_apreds=lapply(dum_zvirus_brts,function(x) x$predict)
zvirus_apreds=do.call(rbind,zvirus_apreds)

## aggregate
zvirus_apreds=data.frame(aggregate(pred~species,data=zvirus_apreds,mean),
                         #aggregate(cpred~species,data=zvirus_apreds,mean)['cpred'], ## holding wos constant
                         aggregate(dum_virus~species,data=virus_apreds,prod)["dum_virus"],
                         aggregate(dum_zvirus~species,data=virus_apreds,prod)["dum_zvirus"])

# ## type
# comp_apreds$type='competence'
# 
# ## apreds
# apreds=rbind.data.frame(pcr_apreds,comp_apreds)

virus_apreds[order(virus_apreds$pred, decreasing=T),]
zvirus_apreds[order(zvirus_apreds$pred, decreasing = T),]

zvirus_apreds$zpred <- zvirus_apreds$pred
zvirus_apreds$pred <- NULL

virus_apreds$vpred <- virus_apreds$pred
virus_apreds$pred <- NULL


# read in data and merge traits
data <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/cleaned dataset 30 cutoff.rds")
apreds <- merge(virus_apreds, zvirus_apreds[c("species","zpred")], by = "species")
apreds <- merge(apreds, data[c("species","Synurbic")], by = "species") # add roost status

write_csv(virus_apreds, "/Volumes/BETKE 2021/bathaus/flat files/virus predictions.csv")
write_csv(zvirus_apreds, "/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus predictions.csv")
write_csv(apreds, "/Volumes/BETKE 2021/bathaus/flat files/all virus predictions.csv")


apreds <- read_csv("/Volumes/BETKE 2021/bathaus/flat files/all virus predictions.csv")
# calculate the number of reservoir species 
# of known reservoirs, how many of them are anthropogenic roosting?
ksyn <- filter(apreds, apreds$dum_virus == 1 & apreds$Synurbic == 1)
known <- filter(apreds, apreds$dum_virus == 1)
nrow(ksyn)/nrow(known)
# 0.6358839

zksyn <- filter(apreds, apreds$dum_zvirus == 1 & apreds$Synurbic == 1)
zknown <- filter(apreds, apreds$dum_zvirus == 1)
nrow(zksyn)/nrow(zknown)
# 0.6541353

# then the number of undetected reservoir species
# of unknown/undetected species, how many are likely to be reservoirs and how many roost in anthropogenic structures?
undet <- filter(apreds, apreds$dum_virus == 0 & vpred > 0.5) %>%
  arrange(desc(vpred)) %>%
  select(!zpred)
syndet <- filter(undet, Synurbic == 1)
nrow(undet) # 37 undetected species, 27 of which are anthropogenic roosting
nrow(syndet) 

# zoonotic 
zundet <- filter(apreds, apreds$dum_zvirus == 0 & zpred > 0.5) %>% 
  arrange(desc(zpred)) %>%
  select(!vpred)
zsyndet <- filter(zundet, Synurbic == 1)
nrow(zundet) # 18 undetected bat species
nrow(zsyndet) # 14 of which are anthropogenic roosting

# save these
write_csv(undet, "/Volumes/BETKE 2021/bathaus/flat files/overall undetected.csv")
write_csv(syndet, "/Volumes/BETKE 2021/bathaus/flat files/overall undetected synurbic.csv")
write_csv(zundet, "/Volumes/BETKE 2021/bathaus/flat files/zoonotic undetected.csv")
write_csv(zsyndet, "/Volumes/BETKE 2021/bathaus/flat files/zoonotic undetected synurbic.csv")

vnames <- undet$species
f <- filter(data, species %in% vnames)
filter(`synurbic and traits only`, species %in% vnames) %>% select(fam) %>% table()

znames <- zundet$species
zf <- filter(data, species %in% znames)
filter(`synurbic and traits only`, species %in% znames) %>% select(fam) %>% table()


### old code I might want to look at later
#### with bat brts
fvirus_fig <- vinfPlot(fvirus_brts, fviurs_df, fvirus_fig, "grey")
fzvirus_fig <- vinfPlot(fzvirus_brts, fzviurs_df, fzvirus_fig, "palegreen3")

# make multi
fvirus <- fvirus_fig[[3]] + 
  labs(x = " ", y = "Relative Importance") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Overall Virus")

fzvirus <- fzvirus_fig[[3]] + 
  labs(x = " ", y = " ") +
  #theme(axis.title.x = element_text(size = 14)) +
  ggtitle("Zoonotic Virus")

# Create the patchwork, dropping the y-axis labels from the plots, and setting
# the margins, this adds the common label
h_patch <- fvirus + fzvirus & ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

# Use the tag label as a y-axis label
png("fvirus_variableinf.png",width=7,height=6,units="in",res=600)
wrap_elements(h_patch) +
  labs(tag = "Relative Importance") +
  theme(
    plot.tag = element_text(size = 10, angle = 90),
    plot.tag.position = "left"
  )
dev.off()

#setwd("/Users/brianabetke/Desktop/Bats and Viruses/ESA 2022/figs")
png("fvirus_variableinf.png",width=7,height=6,units="in",res=600)
fvirus/fzvirus
dev.off()

#### With full bat dataset
virus_fig <- vinfPlot(virus_brts, viurs_df, virus_fig, "dimgrey")
zvirus_fig <- vinfPlot(zvirus_brts, zviurs_df, zvirus_fig, "palegreen3")

# make multi
virus <- fvirus_fig[[3]] + 
  labs(x = " ", y = "Sqrt Relative Importance") +
  theme(axis.title.y = element_text(size = 14, hjust = -14, vjust = 5)) +
  ggtitle("Overall Virus", subtitle = paste("Test AUC =", 
                                            paste(format(round(fvirus_fig[[2]]$avg, 2), nsmall = 2),
                                                  paste(format(round(fvirus_fig[[2]]$se, 2), nsmall = 2),
                                                        sep = " "), sep = " +- ")))
zvirus <- zvirus_fig[[3]] + 
  labs(x = " ", y = " ") +
  theme(axis.title.x = element_text(size = 14)) +
  ggtitle("Zoonotic Virus", subtitle = paste("Test AUC =", 
                                             paste(round(zvirus_fig[[2]]$avg, 2),
                                                   paste(format(round(zvirus_fig[[2]]$se, 2), nsmall = 2),
                                                         sep = " "), sep = " +- ")))

#setwd("/Users/brianabetke/Desktop/Bats and Viruses/ESA 2022/figs")
png("virus_variableinf.png",width=7,height=7,units="in",res=600)
virus/zvirus
dev.off()

# Square root of y
sqrt_virus <- virus + scale_y_sqrt()
sqrt_zvirus <- zvirus + scale_y_sqrt()

png("sqrt_virus_variableinf.png",width=15,height=10,units="in",res=600)
sqrt_virus/sqrt_zvirus
dev.off()





pdp_agg=function(mod, feature){
  
  pdep=plot.gbm(mod[["mod"]], 
                i.var = feature, 
                return.grid = TRUE)
  
  pdep$seed=unique(mod[["seed"]])
  
  pdep$predictor = pdep[feature][,1]
  
  pdep$rank=1:nrow(pdep)
  
  pdep$yhat=pdep$y
  
  return(pdep)
  
}

# aggregate
agg = do.call(rbind, lapply(zoo_prop_brts,function(x) pdp_agg(x, "Synurbic")))
zagg = do.call(rbind, lapply(vrichness_brts,function(x) pdp_agg(x, "Synurbic")))
## get element-wise means
y=with(agg,tapply(yhat,predictor,mean))

## save as mean
#pmean=data.frame(predictor=x,yhat=y)
pmean=data.frame(y)
names(pmean)="yhat"
pmean$predictor=rownames(pmean)
rownames(pmean)=NULL

## make temp data
temp=data
temp$predictor=temp[feature][,1]

## do nothing
agg=agg
pmean=pmean
temp=temp

## get yrange
yrange=range(agg$yhat,pmean$yhat,na.rm=T)

# pull counts for color
df_cat <- as.data.frame(table(temp$predictor))

# fix y axis point
df_cat$ymin <- yrange[1]-0.01

## fix temp to yrange
#temp$yhat=ifelse(temp$predictor==1,max(yrange),min(yrange))

## ggplot with rug
set.seed(1)
ggplot(agg,aes(predictor,yhat,group=seed)) +
  
  ## add individual BRTs
  geom_jitter(size=1,alpha=0.25,colour=pcolor,width=0.1) +
  
  ## add mean
  geom_point(data=pmean,size=2,inherit.aes=F,shape=15,
             aes(predictor,yhat)) +
  
  # # # ## add rug
  # geom_rug(data=df_cat,inherit.aes=F,
  #          aes(Var1, ymin),
  #          sides="b",position="jitter",
  #          colour="grey",alpha=0.25,
  #          na.rm=T)+
  
  geom_point(data=df_cat, inherit.aes = F, shape=23, 
             aes(x=Var1,y=ymin,fill=Freq)) +
  
  ## theme
  theme_bw() +
  theme(axis.text=element_text(size=6),
        axis.title=element_text(size=7)) +
  theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0))) +
  theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0))) +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  labs(x=feature,y="marginal effect") +
  scale_y_continuous(limits=c(yrange[1]-0.01,yrange[2]+0.01),
                     labels=scales::number_format(accuracy=0.01)) +
  scale_fill_continuous(high = "#525252", low = "#D9D9D9", guide="none")
