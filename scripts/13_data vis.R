# 13 - brt visualizations
# babetke@utexas.edu

# clean environment
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

# read in trait data
setwd("/Users/brianabetke/Desktop/bathaus")
data <- readRDS("flat files/log cleaned dataset 30 cutoff.rds")

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
# citation models
cites_brts <- readRDS("flat files/citation brts.rds")
vcites_brts <- readRDS("flat files/virus citation brts.rds")

##### Summary stats of model performance
pull_stats <- function(mod){
  
  if(mod[[1]][["mod"]][["distribution"]][["name"]] == "bernoulli"){
    
    # pull performance stats
    auc <- mean(sapply(mod,function(x) x$testAUC))
    sen <- mean(sapply(mod,function(x) x$sen))
    spec <- mean(sapply(mod,function(x) x$spec))
    se.auc <- std.error(sapply(mod,function(x) x$testAUC)) 
    se.sen <- std.error(sapply(mod,function(x) x$sen))
    se.spec <- std.error(sapply(mod,function(x) x$spec)) 
    
    # df of stats
    df <- data.frame(stat = character(), auc = numeric(), sen = numeric(), spec = numeric())
    df[1,] <- c("mean", auc, sen, spec)
    df[2,] <- c("se", se.auc, se.sen, se.spec)
    return(df)
    
  }else{
    
    # pseudo r2
    test.r2 = mean(sapply(mod,function(x) x$testr2)) 
    se.r2 = std.error(sapply(mod,function(x) x$testr2))
    
    # RMSE
    test.rmse = mean(sapply(mod,function(x) x$testRMSE)) 
    se.rmse = std.error(sapply(mod,function(x) x$testRMSE))
    
    #return
    data.frame(test.r2 = test.r2, se.r2 = se.r2, 
               test.rmse = test.rmse, se.rmse = se.rmse)
  }
}

## Richness models - pseudo R2
pull_stats(vfams_brts)
pull_stats(no_vfams_brts)

## Zoonotic proportion - pseudo R2
pull_stats(zfams_brts)
pull_stats(no_zfams_brts)

## Virus host - AUC
pull_stats(vbinary_brts)
pull_stats(no_vbinary_brts)

## Zoonotic virus host - AUC
pull_stats(zbinary_brts)
pull_stats(no_zbinary_brts)

## Summary stats for citation models
pull_stats(cites_brts)
pull_stats(vcites_brts)

### ttests
## function for extracting data, perform unpaired t test, Cohen's d
# should probably take the with/without labels and model types
# adding data argument to allow for different responses
# Need to allow for addressing negative pseudo R2 values

# break up into 2 functions one where you pull the performance metrics and another 
# that now takes the csv file and makes the figures 
perf_agg=function(mod_with, mod_without, measure){
  
  ## format data
  n=length(sapply(mod_with,function(x) x$measure))
  adata=data.frame(y=c(sapply(mod_with,function(x) x[measure][[1]]),
                       sapply(mod_without,function(x) x[measure][[1]])),
                   response=c(rep('mod_with',n),rep('mod_without',n)),
                   seed=c(sapply(mod_with,function(x) x$seed),
                          sapply(mod_without,function(x) x$seed)))
  rm(n)
  
  # conditionally cutoff negative pseudo R2 values
  if(mod_with[[1]][["mod"]][["distribution"]][["name"]] != "bernoulli") {
    
    # change negatives to 0
    adata <- adata %>% mutate(y = ifelse(y <= 0, 0, y))
    
  }else{
    
    adata <- adata
    
  }
  
  }

# virus models
virus_perf <- perf_agg(vfams_brts, no_vfams_brts, "testr2")
virus_perf_rmse <- perf_agg(vfams_brts, no_vfams_brts, "testRMSE")

# zoonotic models
zfam_perf <- perf_agg(zfams_brts , no_zfams_brts, "testr2")
zfam_perf_rmse <- perf_agg(zfams_brts, no_zfams_brts, "testRMSE")

# virus host models
vbin_perf <- perf_agg(vbinary_brts, no_vbinary_brts, "testAUC")
vbin_SEN <- perf_agg(vbinary_brts, no_vbinary_brts, "sen")
vbin_Spec <- perf_agg(vbinary_brts, no_vbinary_brts, "spec")

# zoonotic host models
zbin_perf <- perf_agg(zbinary_brts, no_zbinary_brts, "testAUC")
zbin_SEN <- perf_agg(zbinary_brts, no_zbinary_brts, "sen")
zbin_Spec <- perf_agg(zbinary_brts, no_zbinary_brts, "spec")

# create t test function
tfun=function(perf_df, fcol){
  
  adata <- perf_df
  
  ## factor
  adata$response=factor(adata$response,levels=c('mod_with','mod_without'))
  
  ## make jitter position
  adata$x=as.numeric(factor(adata$response))
  set.seed(1)
  adata$xj=jitter(adata$x,0.5)
  
  ## fix response (add ifelse statements for synurbic pdps)
  adata$response2=recode(adata$response,
                         "mod_with"="with",
                         "mod_without"="without")
  
  ## t-test
  tsum=t.test(y~response,data=adata,
              alternative='two.sided',
              var.equal=FALSE)  
  
  ## effect size
  csum=cohens_d(y~response,data=adata,paired=F,var.equal=F)

  ## Plot function
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
vrichdata <- tfun(virus_perf, "#E78AC3")

# view stats
vrichdata$tsum
vrichdata$csum

# RMSE
vrichdataRMSE <- tfun(virus_perf_rmse, "#E78AC3")
vrichdataRMSE$tsum
vrichdataRMSE$csum

## Zoo richness
zfamdata <- tfun(zfam_perf, "#FC8D62")

# view stats
zfamdata$tsum
zfamdata$csum

# RMSE
zfamdataRMSE <- tfun(zfam_perf_rmse, "#FC8D62")
zfamdataRMSE$tsum
zfamdataRMSE$csum

## virus reservoir status
# t stats for all metrics
vbinAUC <- tfun(vbin_perf, "#66C2A5")
vbinSEN <- tfun(vbin_SEN, "#66C2A5")
vbinSpec <- tfun(vbin_Spec, "#66C2A5")

# view stats - AUC
vbinAUC$tsum
vbinAUC$csum

# view stats - Sensitivity
vbinSEN$tsum
vbinSEN$csum

# view stats - Specificity
vbinSpec$tsum
vbinSpec$csum

## zoonotic reservoir status
# t stats for all metrics
zbinAUC <- tfun(zbin_perf, "#8DA0CB")
zbinSEN <- tfun(zbin_SEN, "#8DA0CB")
zbinSpec <- tfun(zbin_Spec, "#8DA0CB")

# view stats - AUC
zbinAUC$tsum
zbinAUC$csum

# view stats - Sensitivity
zbinSEN$tsum
zbinSEN$csum

# view stats - Specificity
zbinSpec$tsum
zbinSpec$csum

# all perf. adjustment
ps=c(vrichdata$tsum$p.value,
     zfamdata$tsum$p.value,
     vbinAUC$tsum$p.value,
     vbinSEN$tsum$p.value,
     vbinSpec$tsum$p.value,
     zbinAUC$tsum$p.value,
     zbinSEN$tsum$p.value,
     zbinSpec$tsum$p.value
)
round(p.adjust(ps,method="BH"),4)

# clean out everything but brts and data
rm(list = ls()[!ls() %in% c("data","vfams_brts","no_vfams_brts","zfams_brts","no_zfams_brts","vbinary_brts",
                            "no_vbinary_brts","zbinary_brts", "no_zbinary_brts")])

################### Variable Importance Plots and rankings
# Pull all the relative importance into a dataframe, get the mean, sd, and variation.
# Then create a plot similar to the one I made for the variants 
vinfPlot <- function(data_name, bar_color){
  
  # pull relative importance
  vinf <- lapply(data_name,function(x) x$rinf)
  
  # bind with rbind
  data_vinf <- do.call(rbind,vinf)
  
  # add rankings 
  data_vinf$rank <- rep(1:63, times = 50)
  
  # tidy output
  df_name <- data_vinf %>%
    group_by(var) %>%
    summarize(avg = mean(rel.inf),
              rse = std.error(rel.inf),
              rvar = var(rel.inf)) %>%
    ungroup() %>%
    arrange(desc(avg))
  
  # # calculate average rankings
  # ranks <- data_vinf %>%
  #   group_by(var) %>%
  #   summarize(avg = mean(rank),
  #             med = median(rank),
  #             rse = std.error(rank),
  #             rvar = var(rank)) %>%
  #   ungroup() %>%
  #   arrange(avg)
  
  # Save the raw ranks for boxplots
  ranks <- data_vinf %>% filter(var == "Synurbic")
  
  #Clean up variable names
  df_name$var <- recode(df_name$var,
                         "log_cites" = "Log Citation Count",
                         "log_vcites" = "Log Virus Citation Count",
                         "log_X26.1_GR_Area_km2" = "Log Geographic Area",
                         "X30.1_AET_Mean_mm" = "Mean Monthly AET",
                         "X26.2_GR_MaxLat_dd" = "Maximum Latitude",
                         "X26.3_GR_MinLat_dd" = "Minimum Latitude",
                         "habitat_breadth_n" = "Habitat Breadth",
                         "litters_per_year_n" = "Litters Per Year",
                         "log_adult_body_length_mm" = "Log Adult Body Length",
                         "X28.2_Temp_Mean_01degC" = "Mean Monthly Temperature",
                         "log_X27.2_HuPopDen_Mean_n.km2" = "Log Mean Human Density",
                         "X28.1_Precip_Mean_mm" = "Mean Monthly Precipitation",
                         "litter_size_n" = "Litter Size",
                         "upper_elevation_m" = "Upper Elevation Limit",
                         "disected_by_mountains" = "Disected by Mountains",
                         "log_adult_forearm_length_mm" = "Log Adult Forearm Length",
                         "altitude_breadth_m" = "Altitude Breadth",
                         "X26.4_GR_MidRangeLat_dd" = "Median Latitudinal Range",
                         "foraging_stratum" = "Foraging stratum",
                         "log_adult_mass_g" = "Log Adult Mass",
                         "X30.2_PET_Mean_mm" = "Mean Monthly PET",
                         "det_vfish" = "Diet Fish",
                         "X26.5_GR_MaxLong_dd" = "Maximum Longitude",
                         "fam_RHINOLOPHIDAE" = "Rhinolophidae",
                         "det_diet_breadth_n" = "Diet Breadth",
                         "X26.6_GR_MinLong_dd" = "Minimum Longitude",
                         "det_vend" = "Diet Vend",
                         "log_X27.1_HuPopDen_Min_n.km2" = "Log Min Human Density",
                         "dphy_vertebrate" = "Diet Vertebrate",
                         "log_lower_elevation_m" = "Log Lower Elevation Limit",
                         "det_nect" = "Diet Nectar",
                         "X27.4_HuPopDen_Change" = "Human Density Change",
                         "log_X27.3_HuPopDen_5p_n.km2" = "Log Human Density 5th Percentile",
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
  
  # barplot
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
  return(list(df_name, ranks, fig_name))
}

# Pick colors
library(RColorBrewer)
brewer.pal(n = 8, name = "Set2")

# run for all models with anthropogenic roosting
vrich_fig <- vinfPlot(vfams_brts, "#E78AC3")
zrich_fig <- vinfPlot(zfams_brts, "#FC8D62")
vbinary_fig <- vinfPlot(vbinary_brts, "#66C2A5")
zbinary_fig <- vinfPlot(zbinary_brts, "#8DA0CB")

# save var inf for table?
# write.csv(vrich_fig[[1]], "/Volumes/BETKE 2021/bathaus/flat files/virus richness var inf.csv")
# write.csv(zprop_fig[[1]], "/Volumes/BETKE 2021/bathaus/flat files/zoonotic prop var inf.csv")
# write.csv(vbinary_fig[[1]], "/Volumes/BETKE 2021/bathaus/flat files/virus host var inf.csv")
# write.csv(zbinary_fig[[1]], "/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus host var inf.csv")

# look at plots
vrich_gg <- vrich_fig[[3]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = "Relative Importance (%)") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Virus Family Richness")

zrich_gg <- zrich_fig[[3]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = " ") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Zoonotic Family Richness")

# Create the patchwork, dropping the y-axis labels from the plots, and setting
# the margins, this adds the common label
h_patch <- vrich_gg + zrich_gg & ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

# Use the tag label as a y-axis label
png("figs/figure S3.png",width=7,height=6,units="in",res=600)
wrap_elements(h_patch) +
  labs(tag = "Relative Importance (%)") +
  theme(
    plot.tag = element_text(size = 10),
    plot.tag.position = "bottom"
  )
dev.off()

# binary models
vbinary_gg <- vbinary_fig[[3]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = "Relative Importance (%)") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Virus Host")

zbinary_gg <- zbinary_fig[[3]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = " ") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Zoonotic Host")

# Create the patchwork, dropping the y-axis labels from the plots, and setting
# the margins, this adds the common label
h_patch <- vbinary_gg + zbinary_gg & ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

# Use the tag label as a y-axis label
png("figs/figure S4.png",width=7,height=6,units="in",res=600)
wrap_elements(h_patch) +
  labs(tag = "Relative Importance (%)") +
  theme(
    plot.tag = element_text(size = 10),
    plot.tag.position = "bottom"
  )
dev.off()

# clear from environment
rm(vrich_gg, zrich_gg, vbinary_gg, zbinary_gg, h_patch)

## Rankings
virus <- vrich_fig[[2]]
virus$type <- "virus richness"

zoop <- zrich_fig[[2]]
zoop$type <- "zoonotic richness"

vb <- vbinary_fig[[2]]
vb$type <- "virus host"

zb <- zbinary_fig[[2]]
zb$type <- "zoonotic host"

# rbind
ranks <- rbind(virus, zoop, vb, zb)

# # save
# write.csv(ranks, "/Volumes/BETKE 2021/bathaus/flat files/anthropogenic roost rankings.csv")

# summary
# ranks %>% group_by(type) %>% summarise(med = median(rank), rge = range(rank))
ranks %>% group_by(type) %>% reframe(med = median(rank), rge = range(rank))

# define factor levels
#ranks <- read.csv("/Volumes/BETKE 2021/bathaus/flat files/anthropogenic roost rankings.csv")
ranks$type <- factor(ranks$type, levels = c("virus richness","zoonotic richness","virus host","zoonotic host"))

# ranking plot
# https://schmidtpaul.github.io/DSFAIR/compactletterdisplay.html
# https://statdoe.com/cld-customisation/

# library(multcomp) # reads in mass package which masks select function from dplyr
library(multcompView)

# pairwaise comparisons
mod <- lm(rank ~ type, data = ranks)

# get (adjusted) weight means per group
means <- emmeans::emmeans(object = mod,
                       specs = "type")

# add letters to each mean 
mod_cld <- multcomp::cld(object = means,
                       adjust = "BH",
                       Letters = letters,
                       alpha = 0.05)

# show output
mod_cld

# pull min and max for positioning
lab <- ranks %>% 
  group_by(type) %>% 
  summarise(med = median(rank), min = min(rank), max = max(rank))

# pull letters and add to label df
lab <- merge(lab, mod_cld[c("type",".group")], by = "type")

# then it needs to be added to the plot
png("figs/figure 2.png",width=5,height=4.5,units="in",res=600)
ggplot(ranks, aes(x=type, y=rank, color = type)) + 
  geom_violin() +
  geom_boxplot(width = 0.20) +
  geom_jitter(size=1,alpha=0.25,width=0.2) +
  coord_cartesian(ylim = c(63, 1)) +
  scale_y_reverse() +
  theme_bw() +
  labs(x = "Response", y = "Predictor Rank") +
  scale_color_manual(breaks = c("virus richness","zoonotic richness","virus host","zoonotic host"), 
                     values = c("#E78AC3", "#FC8D62", "#66C2A5", "#8DA0CB")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 8),
        legend.position = "none") +
  geom_text(data = lab, aes(x = type, y = min, label = .group), size = 3, vjust=-2, hjust = 0.5, color = "black")
dev.off()

# Export pairwise and adjusted pvals for supplemental table
tabs <- emmeans::emmeans(mod, pairwise~type, adjust = "BH")
tabs <- summary(tabs$contrasts)

# save test results for supplement 
saveRDS(tabs, "flat files/rank lm.rds")

# clean out everything but brts and data
rm(list = ls()[!ls() %in% c("data","vfams_brts","no_vfams_brts","zfams_brts","no_zfams_brts","vbinary_brts",
                            "no_vbinary_brts","zbinary_brts", "no_zbinary_brts")])

#################### Partial Dependence plots
# so calculate pdep for every model
# Aggregating pdps
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
make_pdp_cont <- function(model,feature,var_name, pcolor) {
  
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
                 size=1,colour="grey",alpha=0.50)+
    
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
    labs(x=var_name,y="marginal effect")+
    scale_y_continuous(labels=scales::number_format(accuracy=0.01))
  
}

## Function for factor pdp plots
make_pdp_fact <- function(model, feature, var_name, pcolor) {
  
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
    
    # geom_point(data=df_cat, inherit.aes = F, shape=23,
    #            aes(x=Var1,y=ymin,fill=Freq)) +

    ## theme
    theme_bw() +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=7)) +
    theme(axis.title.x=element_text(margin=margin(t=5,r=0,b=0,l=0))) +
    theme(axis.title.y=element_text(margin=margin(t=0,r=5,b=0,l=0))) +
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
    labs(x=var_name,y="marginal effect") +
    scale_y_continuous(limits=c(yrange[1]-0.01,yrange[2]+0.01),
                       labels=scales::number_format(accuracy=0.01)) +
    scale_fill_continuous(high = "#525252", low = "#D9D9D9", guide="none")
  
}

# Synurbic pdps (no points - comment out in pdp function)
# MAKE SURE THAT DATA IS READ IN
vsyn <- make_pdp_fact(vfams_brts, "Synurbic", "Anthropogenic Roost", "#E78AC3") 
zsyn <- make_pdp_fact(zfams_brts, "Synurbic", "Anthropogenic Roost", "#FC8D62")
vbsyn <- make_pdp_fact(vbinary_brts, "Synurbic", "Anthropogenic Roost", "#66C2A5")
zbsyn <- make_pdp_fact(zbinary_brts, "Synurbic", "Anthropogenic Roost", "#8DA0CB")

vsyn <- vsyn + 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("no", "yes")) + 
  #ggtitle("Virus Richness") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))
  
zsyn <- zsyn + 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("no", "yes")) + 
  #ggtitle("Zoonotic Proportion") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

vbsyn <- vbsyn + 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("no", "yes")) + 
  #ggtitle("Virus Host") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

zbsyn <- zbsyn + 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("no", "yes")) + 
  #ggtitle("Zoonotic Host") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

png("figs/figure 3.png",width=6,height=4.5,units="in",res=600)
(vsyn + zsyn) / (vbsyn + zbsyn) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 10))
dev.off()

# clean
rm(vsyn, zsyn, vbsyn, zbsyn)

# Remaining pdps
# starting with richness
vit_vrich <- make_pdp_cont(vfams_brts, "log_vcites","Log Virus Citation Count", "#E78AC3")
cit_vrich <- make_pdp_cont(vfams_brts, "log_cites","Log Citation Count", "#E78AC3")
hp_vrich <- make_pdp_cont(vfams_brts, "log_X27.2_HuPopDen_Mean_n.km2", "Log Mean Human Density", "#E78AC3")
aet_vrich <- make_pdp_cont(vfams_brts, "X30.1_AET_Mean_mm", "Mean Monthly AET", "#E78AC3")
mxlon_vrich <- make_pdp_cont(vfams_brts, "X26.5_GR_MaxLong_dd","Maximum Longitude", "#E78AC3")
pet_vrich <- make_pdp_cont(vfams_brts, "X30.2_PET_Mean_mm", "Mean Monthly PET", "#E78AC3")
temp_vrich <- make_pdp_cont(vfams_brts, "X28.2_Temp_Mean_01degC", "Mean Monthly Temperature", "#E78AC3")
hpp_vrich <- make_pdp_cont(vfams_brts, "log_X27.3_HuPopDen_5p_n.km2", "Log Human Density 5th %ile", "#E78AC3")
mxlat_vrich <- make_pdp_cont(vfams_brts, "X26.2_GR_MaxLat_dd", "Maximum Latitude", "#E78AC3")
gr_vrich <- make_pdp_cont(vfams_brts, "log_X26.1_GR_Area_km2", "Log Geographic Area (km2)", "#E78AC3")

png("figs/figure S5.png", width=4,height=6.5,units="in",res=300)
vit_vrich + cit_vrich + hp_vrich + aet_vrich + mxlon_vrich + pet_vrich + temp_vrich + hpp_vrich + mxlat_vrich + gr_vrich + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

# clean
rm(vit_vrich, cit_vrich, hp_vrich, aet_vrich, mxlon_vrich, pet_vrich, temp_vrich, hpp_vrich, mxlat_vrich, gr_vrich)

# zoonotic family richness
vit_zrich <- make_pdp_cont(zfams_brts, "log_vcites","Log Virus Citation Count", "#FC8D62")
cit_zrich <- make_pdp_cont(zfams_brts, "log_cites","Log Citation Count", "#FC8D62")
hp_zrich <- make_pdp_cont(zfams_brts, "log_X27.2_HuPopDen_Mean_n.km2", "Log Mean Human Density", "#FC8D62")
up_zrich <- make_pdp_cont(zfams_brts, "upper_elevation_m", "Upper Elevation Limit", "#FC8D62")
mnlon_zrich <- make_pdp_cont(zfams_brts, "X26.6_GR_MinLong_dd", "Minimum Longitude", "#FC8D62")
gr_zrich <- make_pdp_cont(zfams_brts, "log_X26.1_GR_Area_km2", "Log Geographic Area (km2)", "#FC8D62")
mxlat_zrich <- make_pdp_cont(zfams_brts, "X26.2_GR_MaxLat_dd", "Maximum Latitude", "#FC8D62")
ls_zrich <- make_pdp_cont(zfams_brts, "litter_size_n", "Litter Size", "#FC8D62")
bl_zrich <- make_pdp_cont(zfams_brts, "log_adult_body_length_mm", "Log Adult Body Length", "#FC8D62")
aet_zrich <- make_pdp_cont(zfams_brts, "X30.1_AET_Mean_mm", "Mean Monthly AET", "#FC8D62")

# save
png("figs/figure S6.png", width=4,height=6.5,units="in",res=300)
vit_zrich + cit_zrich + hp_zrich + up_zrich + mnlon_zrich + gr_zrich + mxlat_zrich + ls_zrich + bl_zrich + aet_zrich + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

# clean
rm(vit_zrich, cit_zrich, hp_zrich, up_zrich, mnlon_zrich, gr_zrich, mxlat_zrich, ls_zrich, bl_zrich, aet_zrich)

# virus host status
vit_vbin <- make_pdp_cont(vbinary_brts, "log_vcites","Log Virus Citation Count", "#66C2A5")
cit_vbin <- make_pdp_cont(vbinary_brts, "log_cites","Log Citation Count", "#66C2A5")
gr_vbin <- make_pdp_cont(vbinary_brts, "log_X26.1_GR_Area_km2", "Log Geographic Area (km2)", "#66C2A5")
at_vbin <- make_pdp_cont(vbinary_brts, "X30.1_AET_Mean_mm", "Mean Monthly AET", "#66C2A5")
ml_vbin <- make_pdp_cont(vbinary_brts, "X26.2_GR_MaxLat_dd","Maximum Latitude", "#66C2A5")
mp_vbin <- make_pdp_cont(vbinary_brts, "X30.2_PET_Mean_mm", "Mean Monthly PET", "#66C2A5")
hp_vbin <- make_pdp_cont(vbinary_brts, "log_X27.2_HuPopDen_Mean_n.km2", "Log Mean Human Density", "#66C2A5")
mt_vbin <- make_pdp_cont(vbinary_brts, "X28.2_Temp_Mean_01degC","Mean Monthly Temperature", "#66C2A5")
mi_vbin <- make_pdp_cont(vbinary_brts, "X26.3_GR_MinLat_dd", "Minimum Latitude", "#66C2A5")
fa_vbin <- make_pdp_cont(vbinary_brts, "log_adult_forearm_length_mm", "Log Adult Forearm Length", "#66C2A5")

# save
png("figs/figure S7.png", width=4,height=6.5,units="in",res=300)
vit_vbin + cit_vbin + gr_vbin + at_vbin + ml_vbin + mp_vbin + hp_vbin + mt_vbin + mi_vbin + fa_vbin + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

rm(vit_vbin, cit_vbin, gr_vbin, at_vbin, ml_vbin, mp_vbin, hp_vbin, mt_vbin, mi_vbin, fa_vbin)

# zoonotic host models
vit_zbin <- make_pdp_cont(zbinary_brts, "log_vcites","Log Virus Citation Count", "#8DA0CB")
cit_zbin <- make_pdp_cont(zbinary_brts, "log_cites","Log Citation Count", "#8DA0CB")
gr_zbin <- make_pdp_cont(zbinary_brts, "log_X26.1_GR_Area_km2", "Log Geographic Area (km2)", "#8DA0CB")
at_zbin <- make_pdp_cont(zbinary_brts, "X30.1_AET_Mean_mm", "Mean Monthly AET", "#8DA0CB")
ml_zbin <- make_pdp_cont(zbinary_brts, "X26.2_GR_MaxLat_dd","Maximum Latitude", "#8DA0CB")
mp_zbin <- make_pdp_cont(zbinary_brts,  "X30.2_PET_Mean_mm", "Mean Monthly PET", "#8DA0CB")
mi_zbin <- make_pdp_cont(zbinary_brts, "X26.3_GR_MinLat_dd", "Minimum Latitude", "#8DA0CB")
mt_zbin <- make_pdp_cont(zbinary_brts, "X28.2_Temp_Mean_01degC","Mean Monthly Temperature", "#8DA0CB")
bl_zbin <- make_pdp_cont(zbinary_brts, "log_adult_body_length_mm", "Log Adult Body Length", "#8DA0CB")
mr_zbin <- make_pdp_cont(zbinary_brts, "X26.4_GR_MidRangeLat_dd", "Median Latitudinal Range", "#8DA0CB")

# save
png("figs/figure S8.png", width=4,height=6.5,units="in",res=300)
vit_zbin + cit_zbin + gr_zbin + at_zbin + ml_zbin + mp_zbin + mi_zbin + mt_zbin + bl_zbin + mr_zbin + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

# clean
rm(vit_zbin, cit_zbin, gr_zbin, at_zbin, ml_zbin, mp_zbin, mi_zbin, mt_zbin, bl_zbin, mr_zbin)

############## model predictions
library(ggpubr)
#### Virus richness predictions
vrich_apreds <- lapply(vfams_brts,function(x) x$predict)
vrich_apreds <- do.call(rbind,vrich_apreds)

## aggregate
vrich_apreds  <- data.frame(aggregate(pred~species,data=vrich_apreds,mean),
                        aggregate(cpred~species,data=vrich_apreds,mean)['cpred'])

# am assuming if you want to keep richness in there, you have to merge it back by species?
rich_with <- merge(vrich_apreds, data[c("species","vfam","Synurbic")], by = "species") # add roost status

# label
rich_with$type <- "with"

# without synurbic
# Virus richness predictions
no_vrich_preds <- lapply(no_vfams_brts,function(x) x$predict)
no_vrich_preds <- do.call(rbind,no_vrich_preds)

## aggregate
no_vrich_preds <- data.frame(aggregate(pred~species,data=no_vrich_preds,mean),
                        aggregate(cpred~species,data=no_vrich_preds,mean)['cpred'])

# merge to data to get synurbic
rich_without <- merge(no_vrich_preds, data[c("species","vfam","Synurbic")], by = "species") # add roost status

# type
rich_without$type <- "without"

# either merge for scatterplots or rbind for facets? Both?
all_rich <- bind_rows(rich_with, rich_without)

# facet by known and unknown
all_rich$known <- ifelse(all_rich$vfam > 0, "known", "unknown")

# scatter plot of values with and without
# pivot wider
rich_preds <- all_rich %>% 
  select(species, cpred, type, Synurbic, vfam, known) %>%
  pivot_wider(names_from = type, values_from = cpred)

ggplot(rich_preds, aes((with)^(1/3))) +
  geom_density(aes(fill = Synurbic, color = Synurbic), alpha = 0.5) +
  facet_wrap(~known,ncol=1,strip.position='top',scales="free_y") +
  #scale_x_sqrt() +
  theme_bw() +
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "predicted family richness", fill = "Anthropogenic Roosting", color = "Anthropogenic Roosting") +
  scale_color_manual(labels = c("no","yes","missing"),
                     values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_fill_manual(labels = c("no","yes","missing"),
                    values = c("#8470ff","#9DD866","#A0B1BA")) -> den_rich

# save csv
write.csv(rich_preds, "flat files/richness predictions.csv")

#### zoonotic proportion predictions
zrich_apreds <- lapply(zfams_brts,function(x) x$predict)
zrich_apreds <- do.call(rbind,zrich_apreds)

# aggregate
zrich_apreds <- data.frame(aggregate(pred~species,data=zrich_apreds,mean),
                            aggregate(cpred~species,data=zrich_apreds,mean)['cpred'])

# # get zoonotic prop 
# data <- data %>% 
#   mutate(zoo_prop = zvirus/virus) %>%
#   mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop))

# merge with data
zrich_with <- merge(zrich_apreds, data[c("species","zfam","Synurbic")], by = "species") # add roost status

# label
zrich_with$type <- "with"

# without synurbic
# zoonotic proportion predictions
no_zrich_preds <- lapply(no_zfams_brts,function(x) x$predict)
no_zrich_preds <- do.call(rbind, no_zrich_preds)

## aggregate
no_zrich_preds <- data.frame(aggregate(pred~species,data=no_zrich_preds,mean),
                             aggregate(cpred~species,data=no_zrich_preds,mean)['cpred'])

# merge with data by species - UPDATE TO ZOONOTIC PROP DATA
zrich_without <- merge(no_zrich_preds, data[c("species","zfam","Synurbic")], by = "species") # add roost status

# label
zrich_without$type <- "without"

# rbind
all_zrich <- bind_rows(zrich_with, zrich_without)

# facet by known and unknown
all_zrich$known <- ifelse(all_zrich$zfam > 0, "known", "unknown")

# scatter plot of values with and without
# pivot wider
zrich_preds <- all_zrich %>% 
  select(species, cpred, type, Synurbic, zfam, known) %>%
  pivot_wider(names_from = type, values_from = cpred)

ggplot(zrich_preds, aes((with)^(1/3))) +
  geom_density(aes(fill = Synurbic, color = Synurbic), alpha = 0.5) +
  facet_wrap(~known,ncol=1,strip.position='top',scales="free_y") +
  #scale_x_sqrt(breaks = seq(0,0.75,by=0.05)) +
  theme_bw() +
  theme(legend.position = "top",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "predicted zoonotic family richness", fill = "Anthropogenic Roosting", color = "Anthropogenic Roosting") +
  scale_color_manual(labels = c("no","yes","missing"),
                     values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_fill_manual(labels = c("no","yes","missing"),
                    values = c("#8470ff","#9DD866","#A0B1BA")) -> den_zrich

# save
write.csv(zrich_preds, "flat files/zoonotic family predictions.csv")

## density plots together for supplement
png("figs/figure 4.png", width=6.5,height=5,units="in",res=300)
den_rich + den_zrich + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'A') & 
  theme(legend.position = "bottom")
dev.off()

# # updated facet grid version
# r <- all_rich %>% filter(type == "with") %>% select(-vfam) %>% mutate(outcome = "virus family richness")
# z <- all_zrich %>% filter(type == "with") %>% select(-zfam) %>% mutate(outcome = "zoonotic family richness")
# 
# # bind
# allfams <- rbind(r,z)
# 
# # 
# ggplot(allfams, aes(cpred^(1/3))) + 
#   geom_density(aes(fill=Synurbic,color=Synurbic), alpha = 0.5) +
#   facet_grid(known~outcome) +
#   #scale_x_sqrt(limits = c(0,1), breaks = c(0, 0.2, 0.4, 0.6, 0.8)) +
#   scale_y_sqrt() +
#   theme_bw() +
#   theme(legend.position = "top", 
#         legend.text = element_text(size = 9), 
#         legend.title = element_text(size = 9),
#         legend.key.size = unit(0.4, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   #xlim(0,1) +
#   labs(x = "predicted family richness", fill = "Anthropogenic Roosting", color = "Anthropogenic Roosting") +
#   scale_color_manual(labels = c("no","yes","missing"), 
#                      values = c("#8470ff","#9DD866","#A0B1BA")) +
#   scale_fill_manual(labels = c("no","yes","missing"), 
#                     values = c("#8470ff","#9DD866","#A0B1BA"))

#### average predictions: Overall virus host
# with
vbin_apreds <- lapply(vbinary_brts,function(x) x$predict)
vbin_apreds <- do.call(rbind,vbin_apreds)

## aggregate
vbin_apreds=data.frame(aggregate(pred~species,data=vbin_apreds,mean),
                        aggregate(cpred~species,data=vbin_apreds,mean)['cpred'],
                        aggregate(dum_virus~species,data=vbin_apreds,prod)["dum_virus"])
                        # aggregate(dum_zvirus~species,data=virus_apreds,prod)["dum_zvirus"])

# Merge to species and anthro data
vbin_with <- merge(vbin_apreds, data[c("species","Synurbic")], by = "species") # add roost status

# label
vbin_with$type <- "with"

# without
no_vbin_apreds <- lapply(no_vbinary_brts,function(x) x$predict)
no_vbin_apreds <- do.call(rbind,no_vbin_apreds)

## aggregate
no_vbin_apreds <- data.frame(aggregate(pred~species,data=no_vbin_apreds,mean),
                        aggregate(cpred~species,data=no_vbin_apreds,mean)['cpred'], # holding citations constant
                        aggregate(dum_virus~species,data=no_vbin_apreds,prod)["dum_virus"])


# merge to data to get synurbic
vbin_without <- merge(no_vbin_apreds, data[c("species","Synurbic")], by = "species") # add roost status

# type
vbin_without$type <- "without"

# bind for facets
all_vbin <- bind_rows(vbin_with, vbin_without)

# facet by known and unknown
all_vbin$known <- ifelse(all_vbin$dum_virus == 1, "known", "unknown")

# pivot wider
vbin_preds <- all_vbin %>% 
  select(species, cpred, type, Synurbic, dum_virus) %>%
  pivot_wider(names_from = type, values_from = cpred)

# save
write.csv(vbin_preds, "flat files/virus host predictions.csv")
  
#### average predictions: Zoonotic hosting
zbin_apreds <- lapply(zbinary_brts,function(x) x$predict)
zbin_apreds <- do.call(rbind,zbin_apreds)

# aggregate
zbin_apreds=data.frame(aggregate(pred~species,data=zbin_apreds,mean),
                         aggregate(cpred~species,data=zbin_apreds,mean)['cpred'], ## holding citation count constant
                         aggregate(dum_zvirus~species,data=zbin_apreds,prod)["dum_zvirus"])


# add roost status
zbin_with <- merge(zbin_apreds, data[c("species","Synurbic")], by = "species")

# label
zbin_with$type <- "with"

## Zoonotic without
no_zbin_apreds <- lapply(no_zbinary_brts,function(x) x$predict)
no_zbin_apreds <- do.call(rbind,no_zbin_apreds)

## aggregate
no_zbin_apreds=data.frame(aggregate(pred~species,data=no_zbin_apreds,mean),
                       aggregate(cpred~species,data=no_zbin_apreds,mean)['cpred'], ## holding citation count constant
                       aggregate(dum_zvirus~species,data=no_zbin_apreds,prod)["dum_zvirus"])


# am assuming if you want to keep richness in there, you have to merge it back by species?
zbin_without <- merge(no_zbin_apreds, data[c("species","Synurbic")], by = "species") # add roost status

# label
zbin_without$type <- "without"

# rowbind
all_zbin <- bind_rows(zbin_with, zbin_without)

# pivot wider
zbin_preds <- all_zbin %>% 
  select(species, cpred, type, Synurbic, dum_zvirus) %>%
  pivot_wider(names_from = type, values_from = cpred)

# save
write.csv(zbin_preds, "flat files/zoonotic virus host predictions.csv")

# facet by known and unknown
all_zbin$known <- ifelse(all_zbin$dum_zvirus == 1, "known", "unknown")

# Faceted figure of both host models 
# we need the data we use for the density plots
all_vbin$outcome <- "Virus Hosting"
all_zbin$outcome <- "Zoonotic Virus Hosting"

# filter to with only 
all_vbin %>% filter(type == "with") %>% select(-dum_virus) -> w_vbin
all_zbin %>% filter(type == "with") %>% select(-dum_zvirus) -> w_zbin

# bind
allbin <- rbind(w_vbin,w_zbin)

# ggplot
png("figs/figure 5.png", width=6,height=4.5,units="in",res=300)
ggplot(allbin, aes(cpred)) + 
  geom_density(aes(fill=Synurbic,color=Synurbic), alpha = 0.5) +
  facet_grid(known~outcome) +
  scale_x_sqrt(limits = c(0,1), breaks = c( 0, 0.2, 0.4, 0.6, 0.8)) +
  scale_y_sqrt(breaks = c(0, 2, 4, 6, 8, 10, 12)) +
  theme_bw() +
  theme(legend.position = "top", 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #xlim(0,1) +
  labs(x = expression(paste("predicted probability of hosting (",italic(P),")")), fill = "Anthropogenic Roosting", color = "Anthropogenic Roosting") +
  scale_color_manual(labels = c("no","yes","missing"), 
                     values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_fill_manual(labels = c("no","yes","missing"), 
                    values = c("#8470ff","#9DD866","#A0B1BA"))
dev.off()

##### scatters together for supplemental figure
vfam_scatter <- ggplot(rich_preds, aes(with, without)) +
  geom_point(color = "#E78AC3") +
  theme_bw() +
  labs(x = "with anthropogenic roosting",y = "without anthropogenic roosting")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# scatter plot
zfam_scatter <- ggplot(zrich_preds, aes(with, without)) +
  geom_point(color = "#FC8D62") +
  theme_bw() +
  labs(x = "with anthropogenic roosting", y = "without anthropogenic roosting")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

vbin_scatter <- ggplot(vbin_preds, aes(with, without)) + 
  geom_point(color = "#66C2A5") +
  # geom_smooth(method="lm", color = "grey") +
  theme_bw() +
  labs(x = expression(paste("(",italic(P),") with anthropogenic roosting")), 
       y = expression(paste("(",italic(P),") without anthropogenic roosting"))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0.13,1) +
  ylim(0.13,1) 

zbin_scatter <- ggplot(zbin_preds, aes(with, without)) + 
  geom_point(color = "#8DA0CB") +
  theme_bw() +
  labs(x = expression(paste("(",italic(P),") with anthropogenic roosting")), 
       y = expression(paste("(",italic(P),") without anthropogenic roosting")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0.13,1) +
  ylim(0.13,1)

# add stats
vfam_scatter <- vfam_scatter + stat_cor(method="spearman", cor.coef.name = "rho")
zfam_scatter <- zfam_scatter + stat_cor(method="spearman", cor.coef.name = "rho")
vbin_scatter <- vbin_scatter + stat_cor(method="spearman", cor.coef.name = "rho") 
zbin_scatter <- zbin_scatter + stat_cor(method="spearman", cor.coef.name = "rho") 

# save
png("figs/figure S9.png", width=6.5,height=6,units="in",res=300)
vfam_scatter + zfam_scatter + vbin_scatter + zbin_scatter +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = 'A')
dev.off()
