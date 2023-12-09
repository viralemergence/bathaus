# brt visualizations
# babetke@utexas.edu

#clean environment
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
#data <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/cleaned dataset 30 cutoff.rds")

# lab comp directory
data <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/cleaned dataset 30 cutoff.rds")

###################### BRT Results
# Read in rds files for model outputs from drive
# richness models
vrichness_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/virus with brts.rds")
no_vrichness_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/virus without brts.rds")
# zoonotic proportion models
zoo_prop_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/zoo_prop with brts.rds")
no_zoo_prop_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/zoo_prop without brts.rds")
# virus reservoir
vbinary_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/dum_virus with brts.rds")
no_vbinary_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/dum_virus without brts.rds")
# zoonotic virus reservoir
zbinary_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/dum_zvirus with brts.rds")
no_zbinary_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/dum_zvirus without brts.rds")
# citation models
cites_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/citation brts.rds")
vcites_brts <- readRDS("/Volumes/BETKE 2021/bathaus/flat files/virus citation brts.rds")

##### Summary stats of model performance
## Richness models - pseudo R2
# with
mean(sapply(vrichness_brts,function(x) x$testr2)) # 0.6367943
std.error(sapply(vrichness_brts,function(x) x$testr2)) # 0.01623836

# without
mean(sapply(no_vrichness_brts,function(x) x$testr2)) # 0.6354426
std.error(sapply(no_vrichness_brts,function(x) x$testr2)) # 0.01640279

## Zoonotic proportion - pseudo R2
# with
mean(sapply(zoo_prop_brts,function(x) x$testr2)) # 0.1495972
std.error(sapply(zoo_prop_brts,function(x) x$testr2)) # 0.005234146

# without
mean(sapply(no_zoo_prop_brts,function(x) x$testr2)) # 0.1495524
std.error(sapply(no_zoo_prop_brts,function(x) x$testr2)) # 0.005232724

## Virus reservoir - AUC
# With
# AUC
mean(sapply(vbinary_brts,function(x) x$testAUC)) # 0.9025166
std.error(sapply(vbinary_brts,function(x) x$testAUC)) # 0.001932486
#sen
mean(sapply(vbinary_brts,function(x) x$sen)) # 0.6623377
std.error(sapply(vbinary_brts,function(x) x$sen)) # 0.005751966
# spec
mean(sapply(vbinary_brts,function(x) x$spec)) # 0.9293889
std.error(sapply(vbinary_brts,function(x) x$spec)) # 0.001905445

# without
# AUC
mean(sapply(no_vbinary_brts,function(x) x$testAUC)) # 0.9024412
std.error(sapply(no_vbinary_brts,function(x) x$testAUC)) # 0.001959942
# sen
mean(sapply(no_vbinary_brts,function(x) x$sen)) # 0.6606494
std.error(sapply(no_vbinary_brts,function(x) x$sen)) # 0.005718265
# spec
mean(sapply(no_vbinary_brts,function(x) x$spec)) # 0.9297778
std.error(sapply(no_vbinary_brts,function(x) x$spec)) # 0.001907448

## Zoonotic virus reservoir - AUC
# With
# AUC
mean(sapply(zbinary_brts,function(x) x$testAUC)) # 0.9138378
std.error(sapply(zbinary_brts,function(x) x$testAUC)) # 0.002015201
# sen
mean(sapply(zbinary_brts,function(x) x$sen)) # 0.6172222
std.error(sapply(zbinary_brts,function(x) x$sen)) # 0.006036973
# spec
mean(sapply(zbinary_brts,function(x) x$spec)) # 0.9600493
std.error(sapply(zbinary_brts,function(x) x$spec)) # 0.001355677

# without
# AUC
mean(sapply(no_zbinary_brts,function(x) x$testAUC)) # 0.9138314
std.error(sapply(no_zbinary_brts,function(x) x$testAUC)) # 0.002017853
#sen
mean(sapply(no_zbinary_brts,function(x) x$sen)) # 0.6168519
std.error(sapply(no_zbinary_brts,function(x) x$sen)) # 0.006011903
#spec
mean(sapply(no_zbinary_brts,function(x) x$spec)) # 0.9601478
std.error(sapply(no_zbinary_brts,function(x) x$spec)) # 0.001348789

## need to look at summary stats for citation models
# cites
mean(sapply(cites_brts,function(x) x$testr2)) # 0.1209731
std.error(sapply(cites_brts,function(x) x$testr2)) # 0.07554686
#vcites
mean(sapply(vcites_brts,function(x) x$testr2)) # -0.08667685
std.error(sapply(vcites_brts,function(x) x$testr2)) # 0.1400668

### ttests
## function for extracting data, perform unpaired t test, Cohen's d
# should probably take the with/without labels and model types
# adding data argument to allow for different responses
# Need to allow for addressing negative pseudo R2 values

# break up into 2 functions one where you pull the preformance metrics and another 
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

# call and save as csv files
# virus models
virus_perf <- perf_agg(vrichness_brts, no_vrichness_brts, "testr2")
write_csv(virus_perf, "/Volumes/BETKE 2021/bathaus/flat files/virus performance.csv")
rm(vrichness_brts, no_vrichness_brts)

# zoonotic models
zoop_perf <- perf_agg(zoo_prop_brts , no_zoo_prop_brts, "testr2")
write_csv(zoop_perf, "/Volumes/BETKE 2021/bathaus/flat files/zoonitic proportion performance.csv")
rm(no_zoo_prop_brts, zoo_prop_brts)

# virus host models
vres_perf <- perf_agg(vbinary_brts, no_vbinary_brts, "testAUC")
vres_SEN <- perf_agg(vbinary_brts, no_vbinary_brts, "sen")
vres_Spec <- perf_agg(vbinary_brts, no_vbinary_brts, "spec")
write_csv(vres_perf, "/Volumes/BETKE 2021/bathaus/flat files/virus reservoir performance.csv")
rm(vbinary_brts, no_vbinary_brts)

# zoonotic host models
zoores_perf <- perf_agg(zbinary_brts, no_zbinary_brts, "testAUC")
zoores_SEN <- perf_agg(zbinary_brts, no_zbinary_brts, "sen")
zoores_Spec <- perf_agg(zbinary_brts, no_zbinary_brts, "spec")
write_csv(zoores_perf, "/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus reservoir performance.csv")
rm(zbinary_brts, no_zbinary_brts)

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
vrichdata <- tfun(virus_perf, "#E78AC3")

# view stats
vrichdata$tsum
vrichdata$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = 0.058567, df = 197.98, p-value = 0.9534
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.04416452  0.04686808
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.6367943                 0.6354426

# # A tibble: 1 × 7
# .y.   group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>         <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without 0.00828   100   100 negligible

# # boxplot with significance?
# v_box <-vrichdata[["plot"]] + 
#             labs(x = "Model Type", y = "Model Performance (Pseudo R2)", title = "Virus Richness") +
#             geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
#             geom_text(data = tibble(x=1.5,y = 0.999), 
#                       aes(x=x,y=y, label = paste("T-test: p = ", round(vrichdata$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
#             geom_text(data = tibble(x=1.5,y = 0.925), 
#                       aes(x=x,y=y, label = paste("Cohen's d = ", round(vrichdata$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
#             theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) 
# v_box

## Zoonotic proportion
zpropdata <- tfun(zoop_perf, "#FC8D62")

# view stats
zpropdata$tsum
zpropdata$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = 0.0060513, df = 198, p-value = 0.9952
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.01455050  0.01464007
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.1495972                 0.1495524 

# # A tibble: 1 × 7
# .y.   group1   group2       effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>          <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without 0.000856   100   100 negligible

# # boxplot
# z_box <- zpropdata[["plot"]] + labs(x = "Model Type", y = "Model Performance (Pseudo R2)") +
#             labs(x = "Model Type", y = NULL, title = "Zoonotic Proportion") +
#             geom_line(data = tibble(x=c(1, 2),y = c(0.29, 0.29)), aes(x=x,y=y), inherit.aes = FALSE) +
#             geom_text(data = tibble(x=1.5,y = 0.3), 
#                       aes(x=x,y=y, label = paste("T-test: p = ", round(zpropdata$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
#             geom_text(data = tibble(x=1.5,y = 0.28), 
#                       aes(x=x,y=y, label = paste("Cohen's d = ", round(zpropdata$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
#             theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
# 
# z_box
# 
# v_box + z_box

## virus reservoir status
# need t stats for all metrics
vresAUC <- tfun(vres_perf, "#66C2A5")
vresSEN <- tfun(vres_SEN, "#66C2A5")
vresSpec <- tfun(vres_Spec, "#66C2A5")

# view stats
vresAUC$tsum
vresAUC$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = 0.027393, df = 197.96, p-value = 0.9782
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.005352452  0.005503245
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.9025166                 0.9024412 

# # A tibble: 1 × 7
# .y.   group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>         <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without 0.00387   100   100 negligible

vresSEN$tsum
vresSEN$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = 0.20816, df = 197.99, p-value = 0.8353
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.01430616  0.01768278
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.6623377                 0.6606494 

# # A tibble: 1 × 7
# .y.   group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>         <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without  0.0294   100   100 negligible

vresSpec$tsum
vresSpec$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = -0.14424, df = 198, p-value = 0.8855
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.005705691  0.004927913
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.9293889                 0.9297778 

# # A tibble: 1 × 7
# .y.   group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>         <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without -0.0204   100   100 negligible

# pvalue adjustment for values
ps=c(vresAUC$tsum$p.value,
     vresSEN$tsum$p.value,
     vresSpec$tsum$p.value)
round(p.adjust(ps,method="BH"),4)
# [1] 0.9782 0.9782 0.9782

# # Plot AUC
# # boxplots
# vb_box  <- vresAUC[["plot"]] + labs(x = "Model Type", y = "Model Performance (Test AUC)", title = "Virus Host") +
#               geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
#               geom_text(data = tibble(x=1.5,y = 0.965), 
#                         aes(x=x,y=y, label = paste("T-test: p = ", round(vresAUC$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
#               geom_text(data = tibble(x=1.5,y = 0.955), 
#                         aes(x=x,y=y, label = paste("Cohen's d = ", round(vresAUC$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
#               theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
# vb_box

rm(no_vbinary_brts, vbinary_brts)

## zoonotic reservoir status
# need t stats for all metrics
zresAUC <- tfun(zoores_perf, "#8DA0CB")
zresSEN <- tfun(zoores_SEN, "#8DA0CB")
zresSpec <- tfun(zoores_Spec, "#8DA0CB")

zresAUC$tsum
zresAUC$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = 0.0022392, df = 198, p-value = 0.9982
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.005617415  0.005630186
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.9138378                 0.9138314 

# # A tibble: 1 × 7
# .y.   group1   group2       effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>          <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without 0.000317   100   100 negligible

zresSEN$tsum
zresSEN$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = 0.043471, df = 198, p-value = 0.9654
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.01643094  0.01717169
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.6172222                 0.6168519 

# # A tibble: 1 × 7
# .y.   group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>         <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without 0.00615   100   100 negligible

zresSpec$tsum
zresSpec$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = -0.051519, df = 197.99, p-value = 0.959
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.003869715  0.003672670
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.9600493                 0.9601478 

# # A tibble: 1 × 7
# .y.   group1   group2       effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>          <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without -0.00729   100   100 negligible

# pvalue adjustment for values
ps=c(zresAUC$tsum$p.value,
     zresSEN$tsum$p.value,
     zresSpec$tsum$p.value)
round(p.adjust(ps,method="BH"),4)
# [1] 0.9982 0.9982 0.9982

# # plot AUC
# # boxplots
# zb_box <-  zresAUC[["plot"]] + labs(x = "Model Type", y = NULL, title = "Zoonotic Virus Host") +
#     geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
#     geom_text(data = tibble(x=1.5,y = 0.965), 
#               aes(x=x,y=y, label = paste("T-test: p = ", round(zresAUC$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
#     geom_text(data = tibble(x=1.5,y = 0.955), 
#               aes(x=x,y=y, label = paste("Cohen's d = ", round(zresAUC$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
#     theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
# 
# zb_box

# patch them together and save
png("/Volumes/BETKE 2021/bathaus/figs/ttest performance boxplots.png",width=8,height=8,units="in",res=600)
(v_box + z_box) / (vb_box + zb_box)
dev.off()

# all perf. adjustment
ps=c(vrichdata$tsum$p.value,
     zpropdata$tsum$p.value,
     vresAUC$tsum$p.value,
     vresSEN$tsum$p.value,
     vresSpec$tsum$p.value,
     zresAUC$tsum$p.value,
     zresSEN$tsum$p.value,
     zresSpec$tsum$p.value
)
round(p.adjust(ps,method="BH"),4)
# 0.9982 0.9982 0.9982 0.9982 0.9982 0.9982 0.9982 0.9982

################### Variable Importance Plots and rankings
# Pull all the relative importance into a dataframe, get the mean, sd, and variation.
# Then create a plot similar to the one I made for the variants 
vinfPlot <- function(data_name, bar_color){
  
  # pull relative importance
  vinf <- lapply(data_name,function(x) x$rinf)
  
  # bind with rbind
  data_vinf <- do.call(rbind,vinf)
  
  # add rankings 
  data_vinf$rank <- rep(1:63, times = 100)
  
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

# Pick colors - try to avoid using the same colors as in Fig 1
library(RColorBrewer)
brewer.pal(n = 8, name = "Set2")

# run for all models with anthropogenic roosting
vrich_fig <- vinfPlot(vrichness_brts, "#E78AC3")
zprop_fig <- vinfPlot(zoo_prop_brts, "#FC8D62")
vbinary_fig <- vinfPlot(vbinary_brts, "#66C2A5")
zbinary_fig <- vinfPlot(zbinary_brts, "#8DA0CB")

# save var inf for table?
write.csv(vrich_fig[[1]], "/Volumes/BETKE 2021/bathaus/flat files/virus richness var inf.csv")
write.csv(zprop_fig[[1]], "/Volumes/BETKE 2021/bathaus/flat files/zoonotic prop var inf.csv")
write.csv(vbinary_fig[[1]], "/Volumes/BETKE 2021/bathaus/flat files/virus host var inf.csv")
write.csv(zbinary_fig[[1]], "/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus host var inf.csv")

# look at plots
vrich_gg <- vrich_fig[[3]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = "Relative Importance") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Virus Richness")

zprop_gg <- zprop_fig[[3]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = " ") +
  theme(axis.title.y = element_text(size = 10, hjust = -8, vjust = 5)) +
  ggtitle("Zoonotic Proportion")

# Create the patchwork, dropping the y-axis labels from the plots, and setting
# the margins, this adds the common label
h_patch <- vrich_gg + zprop_gg & ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

# Use the tag label as a y-axis label
png("/Volumes/BETKE 2021/bathaus/figs/figure S3.png",width=7,height=6,units="in",res=600)
wrap_elements(h_patch) +
  labs(tag = "Relative Importance") +
  theme(
    plot.tag = element_text(size = 10),
    plot.tag.position = "bottom"
  )
dev.off()

# binary models
vbinary_gg <- vbinary_fig[[3]] + 
  scale_y_sqrt() + 
  labs(x = " ", y = "Relative Importance") +
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
png("/Volumes/BETKE 2021/bathaus/figs/figure S4.png",width=7,height=6,units="in",res=600)
wrap_elements(h_patch) +
  labs(tag = "Relative Importance") +
  theme(
    plot.tag = element_text(size = 10),
    plot.tag.position = "bottom"
  )
dev.off()

# clear from environment
rm(vrich_gg, zprop_gg, vbinary_gg, zbinary_gg, h_patch)

## Rankings
virus <- vrich_fig[[2]]
virus$type <- "virus richness"

zoop <- zprop_fig[[2]]
zoop$type <- "zoonotic proportion"

vb <- vbinary_fig[[2]]
vb$type <- "virus host"

zb <- zbinary_fig[[2]]
zb$type <- "zoonotic host"

# rbind
ranks <- rbind(virus, zoop, vb, zb)
# try making lolipop plot for synurbic by filtering out anthropogenic roosting values and ggsegement - maybe match colors to var inf plot model colors?
# anthrank <- filter(ranks, var == "Synurbic") %>% 
#   mutate(group = ifelse(type == "virus richness" | type == "zoonotic proportion", "richness", "host"))

# maybe write csv to put these in paper
write.csv(ranks, "/Volumes/BETKE 2021/bathaus/flat files/anthropogenic roost rankings.csv")

# summary
ranks %>% group_by(type) %>% summarise(med = median(rank), rge = range(rank))
# type                  med   rge
# <chr>               <dbl> <int>
# 1 virus host             28    15
# 2 virus host             28    51
# 3 virus richness         37    26
# 4 virus richness         37    44
# 5 zoonotic host          35    28
# 6 zoonotic host          35    47
# 7 zoonotic proportion    42    28
# 8 zoonotic proportion    42    50

# # geom segment
# png("/Volumes/BETKE 2021/bathaus/figs/synurbic ranks.png",width=6.5,height=5,units="in",res=600)
# ggplot(anthrank, aes(x=type, y=ranks, color = type)) +
#   geom_segment(aes(x=type, xend=type, y=ranks, yend=0)) +
#   geom_point(size=5) +
#   coord_flip() +
#   scale_y_reverse() +
#   theme_bw() +
#   labs(x = "Response", y = "Rank") +
#   scale_color_manual(breaks = c("virus richness","zoonotic proportion","virus host","zoonotic host"), 
#                      values = c("#E78AC3", "#FC8D62", "#66C2A5", "#8DA0CB")) +
#   theme(panel.grid.major=element_blank())
# dev.off()
# 

# define factor levels
#ranks <- read.csv("/Volumes/BETKE 2021/bathaus/flat files/anthropogenic roost rankings.csv")
ranks$type <- factor(ranks$type, levels = c("virus richness","zoonotic proportion","virus host","zoonotic host"))

# ranking plot
set.seed(19846)
ggplot(ranks, aes(x=type, y=rank, color = type)) + 
  geom_violin() +
  geom_boxplot(width = 0.20) +
  geom_jitter(size=1,alpha=0.25,width=0.2) +
  # geom_bar(stat = "identity") +
  # geom_errorbar(aes(ymin = med-rse, ymax = med+rse)) +
  #geom_pointrange(aes(ymin = avg-rse, ymax = avg+rse)) +
  #coord_flip() +
  scale_y_reverse() +
  theme_bw() +
  labs(x = "Response", y = "Feature Rank") +
  scale_color_manual(breaks = c("virus richness","zoonotic proportion","virus host","zoonotic host"), 
                     values = c("#E78AC3", "#FC8D62", "#66C2A5", "#8DA0CB")) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 8),
        legend.position = "none") -> ranks_gg

# ttests for ranks
trank <- ranks %>% t_test(rank ~ type, p.adjust.method = "BH")
trank <- trank %>% add_xy_position(x = "type")

png("/Volumes/BETKE 2021/bathaus/figs/figure 2.png",width=5,height=3.5,units="in",res=600)
ranks_gg + ggpubr::stat_pvalue_manual(trank, label = "p.adj.signif", y.position = c(-2,-4,-6,-8,-10,-12), tip.length = 0.01)
dev.off()

# save ttest 
saveRDS(trank, "/Volumes/BETKE 2021/bathaus/flat files/rank ttests.rds")

# clear
rm(virus, zoop, vb, zb, ranks, ranks_gg, trank, vrich_fig, zprop_fig, vbinary_fig, zbinary_fig)

#################### Partial Dependence plots
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
    
    geom_point(data=df_cat, inherit.aes = F, shape=23,
               aes(x=Var1,y=ymin,fill=Freq)) +

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
vsyn <- make_pdp_fact(vrichness_brts, "Synurbic", "Anthropogenic Roost", "#E78AC3") 
zsyn <- make_pdp_fact(zoo_prop_brts, "Synurbic", "Anthropogenic Roost", "#FC8D62")
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

png("/Volumes/BETKE 2021/bathaus/figs/figure 3.png",width=6,height=4.5,units="in",res=600)
(vsyn + zsyn) / (vbsyn + zbsyn) + plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(size = 10))
dev.off()

# clean
rm(vsyn, zsyn, vbsyn, zbsyn)

# Remaining pdps
# starting with richness
vit_rich <- make_pdp_cont(vrichness_brts, "vcites","Virus Citation Count", "#E78AC3")
cit_rich <- make_pdp_cont(vrichness_brts, "cites","Citation Count", "#E78AC3")
gr_rich <- make_pdp_cont(vrichness_brts, "X26.1_GR_Area_km2", "Geographic Area (km2)", "#E78AC3")
bl_rich <- make_pdp_cont(vrichness_brts, "adult_body_length_mm", "Adult Body Length", "#E78AC3")
mx_rich <- make_pdp_cont(vrichness_brts, "X26.2_GR_MaxLat_dd","Maximum Latitude", "#E78AC3")
minln_rich <- make_pdp_cont(vrichness_brts, "X26.6_GR_MinLong_dd", "Minimum Longitude", "#E78AC3")
up_rich <- make_pdp_cont(vrichness_brts, "upper_elevation_m", "Upper Elevation Limit", "#E78AC3")
ls_rich <- make_pdp_cont(vrichness_brts, "litter_size_n", "Litter Size", "#E78AC3")
fa_rich <- make_pdp_cont(vrichness_brts, "adult_forearm_length_mm", "Adult Forearm Length", "#E78AC3")
hb_rich <- make_pdp_cont(vrichness_brts, "habitat_breadth_n", "Habitat Breadth", "#E78AC3")

png("/Volumes/BETKE 2021/bathaus/figs/figure S5.png", width=4,height=6.5,units="in",res=300)
vit_rich + cit_rich + gr_rich + bl_rich + mx_rich + minln_rich + up_rich + ls_rich + fa_rich + hb_rich + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

# clean
rm(vit_rich, cit_rich, gr_rich, bl_rich, mx_rich, minln_rich, up_rich, ls_rich, fa_rich, hb_rich)

# zoonotic proportion
cit_zoop <- make_pdp_cont(zoo_prop_brts, "cites","Citation Count", "#FC8D62")
gr_zoop <- make_pdp_cont(zoo_prop_brts, "X26.1_GR_Area_km2", "Geographic Area (km2)", "#FC8D62")
ml_zoop <- make_pdp_cont(zoo_prop_brts, "X26.2_GR_MaxLat_dd","Maximum Latitude", "#FC8D62")
fa_zoop <- make_pdp_cont(zoo_prop_brts, "adult_forearm_length_mm", "Adult Forearm Length", "#FC8D62")
hp_zoop <- make_pdp_cont(zoo_prop_brts, "X27.2_HuPopDen_Mean_n.km2", "Mean Human Density", "#FC8D62")
vit_zoop <- make_pdp_cont(zoo_prop_brts, "vcites","Virus Citation Count", "#FC8D62")
mp_zoop <- make_pdp_cont(zoo_prop_brts, "X30.2_PET_Mean_mm", "Mean Monthly PET", "#FC8D62")
am_zoop <- make_pdp_cont(zoo_prop_brts, "adult_mass_g", "Adult Mass", "#FC8D62")
bl_zoop <- make_pdp_cont(zoo_prop_brts, "adult_body_length_mm", "Adult Body Length", "#FC8D62")
minln_zoop <- make_pdp_cont(zoo_prop_brts, "X26.6_GR_MinLong_dd", "Minimum Longitude", "#FC8D62")

# save
png("/Volumes/BETKE 2021/bathaus/figs/figure S6.png", width=4,height=6.5,units="in",res=300)
cit_zoop + gr_zoop + ml_zoop + fa_zoop + hp_zoop + vit_zoop + mp_zoop + am_zoop + bl_zoop + minln_zoop + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

# clean
rm(cit_zoop, gr_zoop, ml_zoop, fa_zoop, hp_zoop, vit_zoop, mp_zoop, am_zoop, bl_zoop, minln_zoop)

# virus host status
cit_vres <- make_pdp_cont(vbinary_brts, "cites","Citation Count", "#66C2A5")
vit_vres <- make_pdp_cont(vbinary_brts, "vcites","Virus Citation Count", "#66C2A5")
gr_vres <- make_pdp_cont(vbinary_brts, "X26.1_GR_Area_km2", "Geographic Area (km2)", "#66C2A5")
at_vres <- make_pdp_cont(vbinary_brts, "X30.1_AET_Mean_mm", "Mean Monthly AET", "#66C2A5")
mp_vres <- make_pdp_cont(vbinary_brts, "X30.2_PET_Mean_mm", "Mean Monthly PET", "#66C2A5")
ml_vres <- make_pdp_cont(vbinary_brts, "X26.2_GR_MaxLat_dd","Maximum Latitude", "#66C2A5")
mi_vres <- make_pdp_cont(vbinary_brts, "X26.3_GR_MinLat_dd", "Minimum Latitude", "#66C2A5")
mt_vres <- make_pdp_cont(vbinary_brts, "X28.2_Temp_Mean_01degC","Mean Monthly Temperature", "#66C2A5")
op_vres <- make_pdp_fact(vbinary_brts, "fam_MINIOPTERIDAE", "Miniopteridae", "#66C2A5")
bl_vres <- make_pdp_cont(vbinary_brts, "adult_body_length_mm", "Adult Body Length", "#66C2A5")

# save
png("/Volumes/BETKE 2021/bathaus/figs/figure S7.png", width=4,height=6.5,units="in",res=300)
cit_vres + vit_vres + gr_vres + at_vres + mp_vres + ml_vres + mi_vres + mt_vres + op_vres + bl_vres + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

rm(cit_vres, vit_vres, gr_vres, at_vres, mp_vres, ml_vres, mi_vres, mt_vres, op_vres, bl_vres)

# zoonotic host models
cit_zres <- make_pdp_cont(zbinary_brts, "cites","Citation Count", "#8DA0CB")
vit_zres <- make_pdp_cont(zbinary_brts, "vcites","Virus Citation Count", "#8DA0CB")
gr_zres <- make_pdp_cont(zbinary_brts, "X26.1_GR_Area_km2", "Geographic Area (km2)", "#8DA0CB")
ml_zres <- make_pdp_cont(zbinary_brts, "X26.2_GR_MaxLat_dd","Maximum Latitude", "#8DA0CB")
at_zres <- make_pdp_cont(zbinary_brts, "X30.1_AET_Mean_mm", "Mean Monthly AET", "#8DA0CB")
milon_zres <- make_pdp_cont(zbinary_brts, "X26.6_GR_MinLong_dd", "Minimum Longitude", "#8DA0CB")
fa_zres <- make_pdp_cont(zbinary_brts, "adult_forearm_length_mm", "Adult Forearm Length", "#8DA0CB")
mp_zres <- make_pdp_cont(zbinary_brts,  "X30.2_PET_Mean_mm", "Mean Monthly PET", "#8DA0CB")
bl_zres <- make_pdp_cont(zbinary_brts, "adult_body_length_mm", "Adult Body Length", "#8DA0CB")
up_zres <- make_pdp_cont(zbinary_brts, "upper_elevation_m", "Upper Elevation Limit", "#8DA0CB")

# save
png("/Volumes/BETKE 2021/bathaus/figs/figure S8.png", width=4,height=6.5,units="in",res=300)
cit_zres + vit_zres + gr_zres + ml_zres + at_zres + milon_zres + fa_zres + mp_zres + bl_zres + up_zres + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

# clean
rm(cit_zres, vit_zres, gr_zres, ml_zres, at_zres, milon_zres, fa_zres, mp_zres, bl_zres, up_zres)

############## model predictions
library(ggpubr)
#### Virus richness predictions
vrich_apreds <- lapply(vrichness_brts,function(x) x$predict)
vrich_apreds <- do.call(rbind,vrich_apreds)

## aggregate
vrich_apreds  <- data.frame(aggregate(pred~species,data=vrich_apreds,mean),
                        aggregate(cpred~species,data=vrich_apreds,mean)['cpred'])

# am assuming if you want to keep richness in there, you have to merge it back by species?
rich_with <- merge(vrich_apreds, data[c("species","virus","Synurbic")], by = "species") # add roost status

# label
rich_with$type <- "with"

# without synurbic
# Virus richness predictions
no_vrich_preds <- lapply(no_vrichness_brts,function(x) x$predict)
no_vrich_preds <- do.call(rbind,no_vrich_preds)

## aggregate
no_vrich_preds <- data.frame(aggregate(pred~species,data=no_vrich_preds,mean),
                        aggregate(cpred~species,data=no_vrich_preds,mean)['cpred'])

# merge to data to get synurbic
rich_without <- merge(no_vrich_preds, data[c("species","virus","Synurbic")], by = "species") # add roost status

# type
rich_without$type <- "without"

# either merge for scatterplots or rbind for facets? Both?
all_rich <- bind_rows(rich_with, rich_without)

# # faceted fig
# rich_facet <- ggplot(all_rich, aes(cpred)) + 
#   geom_density(aes(fill = Synurbic, color = Synurbic)) +
#   facet_wrap(~type,ncol=1,strip.position='top',scales="free_y") +
#   theme_bw() +
#   theme(legend.position = "top", 
#         legend.text = element_text(size = 8), 
#         legend.title = element_text(size = 8),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   labs(x = "predicted richness", fill = "Anthropogenic Roosting", color = "Anthropogenic Roosting") +
#   scale_color_manual(labels = c("no","yes","unknown"), 
#                      values = c("#8470ff","#9DD866","#A0B1BA")) +
#   scale_fill_manual(labels = c("no","yes","unknown"), 
#                      values = c("#8470ff","#9DD866","#A0B1BA"))

# scatter plot of values with and without
# pivot wider
rich_preds <- all_rich %>% 
  select(species, cpred, type, Synurbic, virus) %>%
  pivot_wider(names_from = type, values_from = cpred)

# save csv
write.csv(rich_preds, "/Volumes/BETKE 2021/bathaus/flat files/richness predictions.csv")

# # scatter plot
# rich_scatter <- ggplot(rich_preds, aes(with, without)) +
#   geom_point(color = "#E78AC3") +
#   theme_bw() +
#   labs(x = "predicted richness with",y = "predicted richness without")+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())

# save fig
#png("/Volumes/BETKE 2021/bathaus/figs/virus richness preds.png", width=6.5,height=3.5,units="in",res=300)
#rich_facet + rich_scatter + plot_layout(width = c(45, 55), guides = "collect") & theme(legend.position = "top")
#dev.off()

#### zoonotic proportion predictions
zoop_apreds <- lapply(zoo_prop_brts,function(x) x$predict)
zoop_apreds <- do.call(rbind,zoop_apreds)

# aggregate
zoop_apreds <- data.frame(aggregate(back~species,data=zoop_apreds,mean),
                            aggregate(cback~species,data=zoop_apreds,mean)['cback'])

# get zoonotic prop 
data <- data %>% 
  mutate(zoo_prop = zvirus/virus) %>%
  mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop))

# merge with data
zoop_with <- merge(zoop_apreds, data[c("species","zoo_prop","Synurbic")], by = "species") # add roost status

# label
zoop_with$type <- "with"

# without synurbic
# zoonotic proportion predictions
no_zoop_preds <- lapply(no_zoo_prop_brts,function(x) x$predict)
no_zoop_preds <- do.call(rbind,no_zoop_preds)

## aggregate
no_zoop_preds <- data.frame(aggregate(back~species,data=no_zoop_preds,mean),
                             aggregate(cback~species,data=no_zoop_preds,mean)['cback'])

# merge with data by species - UPDATE TO ZOONOTIC PROP DATA
zoop_without <- merge(no_zoop_preds, data[c("species","zoo_prop","Synurbic")], by = "species") # add roost status

# label
zoop_without$type <- "without"

# rbind
all_zoop <- bind_rows(zoop_with, zoop_without)

# faceted fig
zoop_facet <- ggplot(all_zoop, aes(cback)) + 
  geom_density(aes(fill = Synurbic, color = Synurbic)) +
  facet_wrap(~type,ncol=1,strip.position='top',scales="free_y") +
  theme_bw() +
  theme(legend.position = "top", 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "predicted zoonotic proportion", fill = "Anthropogenic Roosting", color = "Anthropogenic Roosting") +
  scale_color_manual(labels = c("no","yes","unknown"), 
                     values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_fill_manual(labels = c("no","yes","unknown"), 
                    values = c("#8470ff","#9DD866","#A0B1BA"))

# scatter plot of values with and without
# pivot wider
zoop_preds <- all_zoop %>% 
  select(species, cback, type, Synurbic, zoo_prop) %>%
  pivot_wider(names_from = type, values_from = cback)

# save
write.csv(zoop_preds, "/Volumes/BETKE 2021/bathaus/flat files/zoonotic proportion predictions.csv")

# # scatter plot
# zoop_scatter <- ggplot(zoop_preds, aes(with, without)) + 
#   geom_point(color = "#FC8D62") +
#   theme_bw() +
#   labs(x = "predicted zoonotic proportion with", y = "predicted zoonotic proportion without")+
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())
# 
# # ggarrange makes them the same size, no problem
# #png("/Volumes/BETKE 2021/bathaus/figs/zoonotic proportion preds.png", width=7,height=3.5,units="in",res=300)
# ggarrange(zoop_facet, zoop_scatter, labels = c("A","B"))
# #dev.off()
# 
# # maybe density plot of knowns instead?
# zoop_preds$known <- ifelse(zoop_preds$zoo_prop == 0, "unknown", "known")
# 
# all_zoop$known <- ifelse(all_zoop$zoo_prop == 0, "unknown", "known")
# 
# ggplot(all_zoop, aes(cback)) +
#   geom_density(aes(fill=known,color=known)) +
#   facet_wrap(~type,ncol=1,strip.position='top',scales="free_y")
# 
# zoop_known <- ggplot(all_zoop, aes(cback)) + 
#   geom_density(aes(fill=known,color=known)) +
#   facet_wrap(~type,ncol=1,strip.position='top',scales="free_y") +
#   theme_bw() +
#   theme(legend.position = "top", 
#         legend.text = element_text(size = 9), 
#         legend.title = element_text(size = 9),
#         legend.key.size = unit(0.4, "cm"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank()) +
#   labs(x = "predicted zoonotic proportion", fill = "Known Associations", color = "Known Associations")
# 
# png("/Volumes/BETKE 2021/bathaus/figs/zoonotic proportion preds Knowns.png", width=7,height=3.5,units="in",res=300)
# ggarrange(zoop_known, zoop_scatter, labels = c("A","B"))
# dev.off()

#### average predictions: Overall virus reservoir
# with
vres_apreds <- lapply(vbinary_brts,function(x) x$predict)
vres_apreds <- do.call(rbind,vres_apreds)

## aggregate
vres_apreds=data.frame(aggregate(pred~species,data=vres_apreds,mean),
                        aggregate(cpred~species,data=vres_apreds,mean)['cpred'],
                        aggregate(dum_virus~species,data=vres_apreds,prod)["dum_virus"])
                        # aggregate(dum_zvirus~species,data=virus_apreds,prod)["dum_zvirus"])

# am assuming if you want to keep richness in there, you have to merge it back by species?
vres_with <- merge(vres_apreds, data[c("species","Synurbic")], by = "species") # add roost status

# label
vres_with$type <- "with"

# without
no_vres_apreds <- lapply(no_vbinary_brts,function(x) x$predict)
no_vres_apreds <- do.call(rbind,no_vres_apreds)

## aggregate
no_vres_apreds <- data.frame(aggregate(pred~species,data=no_vres_apreds,mean),
                        aggregate(cpred~species,data=no_vres_apreds,mean)['cpred'], # holding citations constant
                        aggregate(dum_virus~species,data=no_vres_apreds,prod)["dum_virus"])


# merge to data to get synurbic
vres_without <- merge(no_vres_apreds, data[c("species","Synurbic")], by = "species") # add roost status

# type
vres_without$type <- "without"

# either merge for scatterplots or rbind for facets? Both?
all_vres <- bind_rows(vres_with, vres_without)

# facet by known and unknown
all_vres$known <- ifelse(all_vres$dum_virus == 1, "unknown", "known")

# plot for with only
all_vres %>% 
  filter(type == "with") %>%
ggplot(aes(cpred)) + 
  geom_density(aes(fill=Synurbic,color=Synurbic), alpha = 0.5) +
  facet_wrap(~known,ncol=1,strip.position='top',scales="free_y") +
  theme_bw() +
  theme(legend.position = "top", 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = expression(paste("predicted probablity of hosting (",italic(P),")")), fill = "Anthropogenic Roosting", color = "Anthropogenic Roosting") +
  scale_color_manual(labels = c("no","yes","unknown"), 
                     values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_fill_manual(labels = c("no","yes","unknown"), 
                    values = c("#8470ff","#9DD866","#A0B1BA")) -> vres_gg

# scatter plot of values with and without
# pivot wider
vres_preds <- all_vres %>% 
  select(species, cpred, type, Synurbic, dum_virus) %>%
  pivot_wider(names_from = type, values_from = cpred)

# save
write.csv(vres_preds, "/Volumes/BETKE 2021/bathaus/flat files/virus host predictions.csv")

# scatter plot
vres_scatter <- ggplot(vres_preds, aes(with, without)) + 
  geom_point(color = "#66C2A5") +
  theme_bw() +
  labs(x = expression(paste("(",italic(P),") with anthropogenic roosting")), 
       y = expression(paste("(",italic(P),") without anthropogenic roosting"))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
# save figs 
png("/Volumes/BETKE 2021/bathaus/figs/figure 4.png", width=7,height=3.5,units="in",res=300)
ggarrange(vres_scatter, vres_gg, labels = c("A","B"))
dev.off()

# cor test for virus hosts
cor.test(vres_preds$with, vres_preds$without, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  vres_preds$with and vres_preds$without
# t = 752.2, df = 1277, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9987429 0.9989905
# sample estimates:
#       cor 
# 0.9988734  

#### average predictions: Zoonotic
zres_apreds <- lapply(zbinary_brts,function(x) x$predict)
zres_apreds <- do.call(rbind,zres_apreds)

# aggregate
zres_apreds=data.frame(aggregate(pred~species,data=zres_apreds,mean),
                         aggregate(cpred~species,data=zres_apreds,mean)['cpred'], ## holding citation count constant
                         aggregate(dum_zvirus~species,data=zres_apreds,prod)["dum_zvirus"])


# add roost status
zres_with <- merge(zres_apreds, data[c("species","Synurbic")], by = "species")

# label
zres_with$type <- "with"

## Zoonotic without
no_zres_apreds <- lapply(no_zbinary_brts,function(x) x$predict)
no_zres_apreds <- do.call(rbind,no_zres_apreds)

## aggregate
no_zres_apreds=data.frame(aggregate(pred~species,data=no_zres_apreds,mean),
                       aggregate(cpred~species,data=no_zres_apreds,mean)['cpred'], ## holding citation count constant
                       aggregate(dum_zvirus~species,data=no_zres_apreds,prod)["dum_zvirus"])


# am assuming if you want to keep richness in there, you have to merge it back by species?
zres_without <- merge(no_zres_apreds, data[c("species","Synurbic")], by = "species") # add roost status

# label
zres_without$type <- "without"

# rowbind
all_zres <- bind_rows(zres_with, zres_without)

# scatter plot of values with and without
# pivot wider
zres_preds <- all_zres %>% 
  select(species, cpred, type, Synurbic, dum_zvirus) %>%
  pivot_wider(names_from = type, values_from = cpred)

# save
write.csv(zres_preds, "/Volumes/BETKE 2021/bathaus/flat files/zoonotic virus host predictions.csv")

# plots
# scatter plot
# facet by known and unknown
all_zres$known <- ifelse(all_zres$dum_zvirus == 1, "unknown", "known")

# plot for with only
all_zres %>% 
  filter(type == "with") %>%
  ggplot(aes(cpred)) + 
  geom_density(aes(fill=Synurbic,color=Synurbic), alpha = 0.5) +
  facet_wrap(~known,ncol=1,strip.position='top',scales="free_y") +
  theme_bw() +
  theme(legend.position = "top", 
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0,1) +
  labs(x = expression(paste("predicted probablity of hosting zoonotic virus (",italic(P),")")), fill = "Anthropogenic Roosting", color = "Anthropogenic Roosting") +
  scale_color_manual(labels = c("no","yes","unknown"), 
                     values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_fill_manual(labels = c("no","yes","unknown"), 
                    values = c("#8470ff","#9DD866","#A0B1BA")) -> zres_gg

zres_scatter <- ggplot(zres_preds, aes(with, without)) + 
  geom_point(color = "#8DA0CB") +
  theme_bw() +
  labs(x = expression(paste("(",italic(P),") with anthropogenic roosting")), 
       y = expression(paste("(",italic(P),") without anthropogenic roosting")))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(0.13,1) +
  ylim(0.13,1)

# save
png("/Volumes/BETKE 2021/bathaus/figs/figure 5.png", width=7.25,height=3.5,units="in",res=300)
ggarrange(zres_scatter, zres_gg, labels = c("A","B"))
dev.off()

# cor test
cor.test(zres_preds$with, zres_preds$without, method = "pearson")
# Pearson's product-moment correlation
# 
# data:  zres_preds$with and zres_preds$without
# t = 4810.7, df = 1277, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9999692 0.9999753
# sample estimates:
#       cor 
# 0.9999724 
