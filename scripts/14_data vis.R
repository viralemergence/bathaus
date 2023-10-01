# data vis - descriptive and brts?

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
data <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/cleaned dataset 30 cutoff.rds")

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
mean(sapply(vrichness_brts,function(x) x$testr2)) # 0.5572603
std.error(sapply(vrichness_brts,function(x) x$testr2)) # 0.03571409
# without
mean(sapply(no_vrichness_brts,function(x) x$testr2)) # 0.5544271
std.error(sapply(no_vrichness_brts,function(x) x$testr2)) # 0.0363453

## Zoonotic proportion - pseudo R2
# with
mean(sapply(zoo_prop_brts,function(x) x$testr2)) #  0.1333257
std.error(sapply(zoo_prop_brts,function(x) x$testr2)) # 0.005387098
# without
mean(sapply(no_zoo_prop_brts,function(x) x$testr2)) # 0.1331381
std.error(sapply(no_zoo_prop_brts,function(x) x$testr2)) # 0.005441455

## Virus reservoir - AUC
# With
# AUC
mean(sapply(vbinary_brts,function(x) x$testAUC)) # 0.9031327
std.error(sapply(vbinary_brts,function(x) x$testAUC)) # 0.001730591
#sen
mean(sapply(vbinary_brts,function(x) x$sen)) # 0.6742105
std.error(sapply(vbinary_brts,function(x) x$sen)) # 0.004808659
# spec
mean(sapply(vbinary_brts,function(x) x$spec)) # 0.9279444
std.error(sapply(vbinary_brts,function(x) x$spec)) # 0.00193009

# without
# AUC
mean(sapply(no_vbinary_brts,function(x) x$testAUC)) # 0.9032654
std.error(sapply(no_vbinary_brts,function(x) x$testAUC)) # 0.001741602
# sen
mean(sapply(no_vbinary_brts,function(x) x$sen)) # 0.6742105
std.error(sapply(no_vbinary_brts,function(x) x$sen)) # 0.004764818
# spec
mean(sapply(no_vbinary_brts,function(x) x$spec)) # 0.9287222
std.error(sapply(no_vbinary_brts,function(x) x$spec)) # 0.001819678

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
mean(sapply(cite_brts,function(x) x$testr2)) # 0.1209731
std.error(sapply(cite_brts,function(x) x$testr2)) # 0.07554686
#vcites
mean(sapply(vcite_brts,function(x) x$testr2)) # -0.08667685
std.error(sapply(vcite_brts,function(x) x$testr2)) # 0.1400668

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
  
  # conditionally cut negative pseudo R2 values to 0
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
# t = 0.071678, df = 197.98, p-value = 0.9429
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.05606260  0.06029181
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.5819087                 0.5797941 

# # A tibble: 1 × 7
#    .y.  group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>        <dbl>   <int> <int> <ord>     
# 1   y   mod_with mod_without  0.0101   100   100 negligible

# boxplot with significance?
v_box <-vrichdata[["plot"]] + 
            labs(x = "Model Type", y = "Model Performance (Pseudo R2)", title = "Virus Richness") +
            geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
            geom_text(data = tibble(x=1.5,y = 0.999), 
                      aes(x=x,y=y, label = paste("T-test: p = ", round(vrichdata$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
            geom_text(data = tibble(x=1.5,y = 0.925), 
                      aes(x=x,y=y, label = paste("Cohen's d = ", round(vrichdata$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
            theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold")) 
v_box

## Zoonotic proportion
zpropdata <- tfun(zoop_perf, "#FC8D62")

# view stats
zpropdata$tsum
zpropdata$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = 0.025077, df = 197.98, p-value = 0.98
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.01481379  0.01519540
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.1334578                 0.1332670 

# # A tibble: 1 × 7
# .y.   group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>         <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without 0.00355   100   100 negligible

# boxplot
z_box <- zpropdata[["plot"]] + labs(x = "Model Type", y = "Model Performance (Pseudo R2)") +
            labs(x = "Model Type", y = NULL, title = "Zoonotic Proportion") +
            geom_line(data = tibble(x=c(1, 2),y = c(0.29, 0.29)), aes(x=x,y=y), inherit.aes = FALSE) +
            geom_text(data = tibble(x=1.5,y = 0.3), 
                      aes(x=x,y=y, label = paste("T-test: p = ", round(zpropdata$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
            geom_text(data = tibble(x=1.5,y = 0.28), 
                      aes(x=x,y=y, label = paste("Cohen's d = ", round(zpropdata$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
            theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

z_box

v_box + z_box

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
# t = -0.054038, df = 197.99, p-value = 0.957
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.004974419  0.004709068
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.9031327                 0.9032654 

# # A tibble: 1 × 7
# .y.   group1   group2       effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>          <dbl> <int> <int> <ord>     
#   1 y     mod_with mod_without -0.00764   100   100 negligible

vresSEN$tsum
vresSEN$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = 0, df = 197.98, p-value = 1
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.01334966  0.01334966
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.6742105                 0.6742105 

# A tibble: 1 × 7
# .y.   group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>         <dbl> <int> <int> <ord>     
# 1 y     mod_with mod_without       0   100   100 negligible

vresSpec$tsum
vresSpec$csum

# Welch Two Sample t-test
# 
# data:  y by response
# t = -0.29321, df = 197.32, p-value = 0.7697
# alternative hypothesis: true difference in means between group mod_with and group mod_without is not equal to 0
# 95 percent confidence interval:
#   -0.006008933  0.004453377
# sample estimates:
#   mean in group mod_with mean in group mod_without 
# 0.9279444                 0.9287222

# # A tibble: 1 × 7
# .y.   group1   group2      effsize    n1    n2 magnitude 
# * <chr> <chr>    <chr>         <dbl> <int> <int> <ord>     
# 1 y     mod_with mod_without -0.0415   100   100 negligible

# pvalue adjustment for values
ps=c(vresAUC$tsum$p.value,
     vresSEN$tsum$p.value,
     vresSpec$tsum$p.value)
round(p.adjust(ps,method="BH"),4)
# [1] 1 1 1

# Plot AUC
# boxplots
vb_box  <- vresAUC[["plot"]] + labs(x = "Model Type", y = "Model Performance (Test AUC)", title = "Virus Host") +
              geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
              geom_text(data = tibble(x=1.5,y = 0.965), 
                        aes(x=x,y=y, label = paste("T-test: p = ", round(vresAUC$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
              geom_text(data = tibble(x=1.5,y = 0.955), 
                        aes(x=x,y=y, label = paste("Cohen's d = ", round(vresAUC$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
              theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
vb_box

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

# plot AUC
# boxplots
zb_box <-  zresAUC[["plot"]] + labs(x = "Model Type", y = NULL, title = "Zoonotic Virus Host") +
    geom_line(data = tibble(x=c(1, 2),y = c(0.96, 0.96)), aes(x=x,y=y), inherit.aes = FALSE) +
    geom_text(data = tibble(x=1.5,y = 0.965), 
              aes(x=x,y=y, label = paste("T-test: p = ", round(zresAUC$tsum$p.value, 4), sep = "")), inherit.aes = FALSE) +
    geom_text(data = tibble(x=1.5,y = 0.955), 
              aes(x=x,y=y, label = paste("Cohen's d = ", round(zresAUC$csum$effsize, 4), sep = "")), inherit.aes = FALSE) +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))

zb_box

# patch them together and save
png("/Volumes/BETKE 2021/bathaus/figs/ttest performance boxplots.png",width=8,height=8,units="in",res=600)
(v_box + z_box) / (vb_box + zb_box)
dev.off()

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
  
  # calculate average rankings
  ranks <- data_vinf %>%
    group_by(var) %>%
    summarize(avg = mean(rank),
              rse = std.error(rank),
              rvar = var(rank)) %>%
    ungroup() %>%
    arrange(avg)
  
  # # pull relative importance
  # vinf <- lapply(data_name,function(x) x$rinf)
  # 
  # # bind with rbind
  # data_vinf <- do.call(rbind,vinf)
  # 
  # # tidy output
  # df_name <- data_vinf %>%
  #   group_by(var) %>%
  #   summarize(avg = mean(rel.inf),
  #             rse = std.error(rel.inf),
  #             rvar = var(rel.inf)) %>%
  #   ungroup() %>%
  #   arrange(desc(avg))
  
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
png("/Volumes/BETKE 2021/bathaus/figs/richnes_variableinf.png",width=7,height=6,units="in",res=600)
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
virus <- vrich_fig[[2]]
#virus$ranks <- 1:nrow(virus)
virus$type <- "virus richness"

zoop <- zprop_fig[[2]]
#zoop$ranks <- 1:nrow(zoop)
zoop$type <- "zoonotic proportion"

vb <- vbinary_fig[[2]]
#vb$ranks <- 1:nrow(vb)
vb$type <- "virus host"

zb <- zbinary_fig[[2]]
#zb$ranks <- 1:nrow(zb)
zb$type <- "zoonotic host"

# merge into one dataset by variable
# try making lolipop plot for synurbic by filtering out anthropogenic roosting values and ggsegement - maybe match colors to var inf plot model colors?
ranks <- rbind(virus, zoop, vb, zb)
anthrank <- filter(ranks, var == "Synurbic") %>% 
  mutate(group = ifelse(type == "virus richness" | type == "zoonotic proportion", "richness", "host"))

# maybe write csv to put these in paper
write.csv(anthrank, "/Volumes/BETKE 2021/bathaus/flat files/anthropogenic roost rankings.csv")

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
anthrank$type <- factor(anthrank$type, levels = c("virus richness","zoonotic proportion","virus host","zoonotic host"))

# alternative plot - reverse barplot with error bars
png("/Volumes/BETKE 2021/bathaus/figs/synurbic ranks.png",width=4.5,height=3,units="in",res=600)
ggplot(anthrank, aes(x=type, y=avg, fill = type)) + 
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg-rse, ymax = avg+rse)) +
  #geom_pointrange(aes(ymin = avg-rse, ymax = avg+rse)) +
  #coord_flip() +
  scale_y_reverse() +
  theme_bw() +
  labs(x = "Response", y = "Rank") +
  scale_fill_manual(breaks = c("virus richness","zoonotic proportion","virus host","zoonotic host"), 
                     values = c("#E78AC3", "#FC8D62", "#66C2A5", "#8DA0CB")) +
  theme(# axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
    legend.position = "none") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 8))
dev.off()

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

# Synurbic pdps
# MAKE SURE THAT DATA IS READ IN
vsyn <- make_pdp_fact(vrichness_brts, "Synurbic", "Anthropogenic Roost", "#E78AC3") 
zsyn <- make_pdp_fact(zoo_prop_brts, "Synurbic", "Anthropogenic Roost", "#FC8D62")
vbsyn <- make_pdp_fact(vbinary_brts, "Synurbic", "Anthropogenic Roost", "#66C2A5")
zbsyn <- make_pdp_fact(zbinary_brts, "Synurbic", "Anthropogenic Roost", "#8DA0CB")

vsyn <- vsyn + 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("no", "yes")) + 
  ggtitle("Virus Richness") + 
  theme(plot.title = element_text(hjust = 0.5, size = 10))
  
zsyn <- zsyn + 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("no", "yes")) + 
  ggtitle("Zoonotic Proportion") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

vbsyn <- vbsyn + 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("no", "yes")) + 
  ggtitle("Virus Host") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

zbsyn <- zbsyn + 
  scale_x_discrete(breaks=c("0","1"),
                   labels=c("no", "yes")) + 
  ggtitle("Zoonotic Host") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

png("/Volumes/BETKE 2021/bathaus/figs/synurbic pdps.png",width=6,height=4.5,units="in",res=600)
(vsyn + zsyn) / (vbsyn + zbsyn)
dev.off()

# ttest and cohens D for synurbic pdps
pdp_tfun <- function(model){
  
  agg = do.call(rbind, lapply(model,function(x) pdp_agg(x,"Synurbic")))
  
  ## t-test
  tsum=t.test(y~Synurbic,data=agg,
              alternative='two.sided',
              var.equal=F,paired=F)  
  
  ## effect size
  csum=cohens_d(y~Synurbic,data=agg,paired=F,var.equal=F)
  
  return(list(tsum, csum))
  
}

pdp_tfun(vrichness_brts)
# [[1]]
# 
# Welch Two Sample t-test
# 
# data:  y by Synurbic
# t = 0.25401, df = 198, p-value = 0.7997
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#   -0.05079995  0.06582183
# sample estimates:
#   mean in group 0 mean in group 1 
# 3.405688        3.398177 
# 
# 
# [[2]]
# # A tibble: 1 × 7
# .y.   group1 group2 effsize    n1    n2 magnitude 
# * <chr> <chr>  <chr>    <dbl> <int> <int> <ord>     
#   1 y     0      1       0.0359   100   100 negligible

pdp_tfun(zoo_prop_brts)
# [[1]]
# 
# Welch Two Sample t-test
# 
# data:  y by Synurbic
# t = -1.6757, df = 197.99, p-value = 0.09538
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#   -0.0029677012  0.0002411034
# sample estimates:
#   mean in group 0 mean in group 1 
# 0.1341661       0.1355294 

# [[2]]
# # A tibble: 1 × 7
# .y.   group1 group2 effsize    n1    n2 magnitude
# * <chr> <chr>  <chr>    <dbl> <int> <int> <ord>    
#   1 y     0      1       -0.237   100   100 small   

pdp_tfun(vbinary_brts)

# [[1]]
# 
# Welch Two Sample t-test
# 
# data:  y by Synurbic
# t = -7.8138, df = 197.98, p-value = 3.194e-13
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#   -0.04252711 -0.02538708
# sample estimates:
#   mean in group 0 mean in group 1 
# 0.2394832       0.2734403 
# 
# 
# [[2]]
# # A tibble: 1 × 7
# .y.   group1 group2 effsize    n1    n2 magnitude
# * <chr> <chr>  <chr>    <dbl> <int> <int> <ord>    
#   1 y     0      1        -1.11   100   100 large    

pdp_tfun(zbinary_brts)

# [[1]]
# 
# Welch Two Sample t-test
# 
# data:  y by Synurbic
# t = -1.7076, df = 197.93, p-value = 0.08927
# alternative hypothesis: true difference in means between group 0 and group 1 is not equal to 0
# 95 percent confidence interval:
#   -0.0043504965  0.0003125891
# sample estimates:
#   mean in group 0 mean in group 1 
# 0.1533488       0.1553677 
# 
# 
# [[2]]
# # A tibble: 1 × 7
# .y.   group1 group2 effsize    n1    n2 magnitude
# * <chr> <chr>  <chr>    <dbl> <int> <int> <ord>    
#   1 y     0      1       -0.241   100   100 small    

# Remaining pdps
# starting with richness
cit_rich <- make_pdp_cont(vrichness_brts, "cites","Citation Count", "#E78AC3")
vit_rich <- make_pdp_cont(vrichness_brts, "vcites","Virus Citation Count", "#E78AC3")
gr_rich <- make_pdp_cont(vrichness_brts, "X26.1_GR_Area_km2", "Geographic Area (km2)", "#E78AC3")
mx_rich <- make_pdp_cont(vrichness_brts, "X26.2_GR_MaxLat_dd","Maximum Latitude", "#E78AC3")
hp_rich <- make_pdp_cont(vrichness_brts, "X27.2_HuPopDen_Mean_n.km2", "Mean Human Density", "#E78AC3")
hb_rich <- make_pdp_cont(vrichness_brts, "habitat_breadth_n", "Habitat Breadth", "#E78AC3")
ml_rich <- make_pdp_cont(vrichness_brts, "X26.3_GR_MinLat_dd", "Minimum Latitude", "#E78AC3")
mp_rich <- make_pdp_cont(vrichness_brts, "X30.2_PET_Mean_mm", "Mean Monthly PET", "#E78AC3")
ls_rich <- make_pdp_cont(vrichness_brts, "litter_size_n", "Litter Size", "#E78AC3")
mlon_rich <- make_pdp_cont(vrichness_brts, "X26.5_GR_MaxLong_dd", "Maximum Longitude", "#E78AC3")

png("/Volumes/BETKE 2021/bathaus/figs/richness pdps.png", width=4,height=6.5,units="in",res=300)
vit_rich + cit_rich + gr_rich + mx_rich + hp_rich + hb_rich + ml_rich + mp_rich + ls_rich + mlon_rich + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

# zoonotic proportion
cit_zoop <- make_pdp_cont(zoo_prop_brts, "cites","Citation Count", "#FC8D62")
vit_zoop <- make_pdp_cont(zoo_prop_brts, "vcites","Virus Citation Count", "#FC8D62")
ml_zoop <- make_pdp_cont(zoo_prop_brts, "X26.2_GR_MaxLat_dd","Maximum Latitude", "#FC8D62")
fa_zoop <- make_pdp_cont(zoo_prop_brts, "adult_forearm_length_mm", "Adult Forearm Length", "#FC8D62")
hp_zoop <- make_pdp_cont(zoo_prop_brts, "X27.2_HuPopDen_Mean_n.km2", "Mean Human Density", "#FC8D62")
mp_zoop <- make_pdp_cont(zoo_prop_brts, "X30.2_PET_Mean_mm", "Mean Monthly PET", "#FC8D62")
ls_zoop <- make_pdp_cont(zoo_prop_brts, "litter_size_n", "Litter Size", "#FC8D62")
bl_zoop <- make_pdp_cont(zoo_prop_brts, "adult_body_length_mm", "Adult Body Length", "#FC8D62")
mt_zoop <- make_pdp_cont(zoo_prop_brts, "X28.2_Temp_Mean_01degC","Mean Monthly Temperature", "#FC8D62")
nt_zoop <- make_pdp_fact(zoo_prop_brts, "Neotropical", "Neotropical", "#FC8D62") 

# save
png("/Volumes/BETKE 2021/bathaus/figs/zoonotic proportion pdps.png", width=4,height=6.5,units="in",res=300)
cit_zoop + vit_zoop + ml_zoop + fa_zoop + hp_zoop + mp_zoop + ls_zoop + bl_zoop + mt_zoop + nt_zoop + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

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

#save
png("/Volumes/BETKE 2021/bathaus/figs/virus host pdps.png", width=4,height=6.5,units="in",res=300)
cit_vres + vit_vres + gr_vres + at_vres + mp_vres + ml_vres + mi_vres + mt_vres + op_vres + bl_vres + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

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
png("/Volumes/BETKE 2021/bathaus/figs/zoonotic virus host pdps.png", width=4,height=6.5,units="in",res=300)
cit_zres + vit_zres + gr_zres + ml_zres + at_zres + milon_zres + fa_zres + mp_zres + bl_zres + up_zres + plot_layout(nrow = 5, ncol = 2, byrow = TRUE)
dev.off()

# one where the top 5 for all of them in the same plot?
png("/Volumes/BETKE 2021/bathaus/figs/All pdps.png", width=8,height=5,units="in",res=300)
vit_rich + cit_rich + gr_rich + mx_rich + hp_rich +
cit_zoop + vit_zoop + ml_zoop + fa_zoop + hp_zoop +
cit_vres + vit_vres + gr_vres + at_vres + mp_vres +
cit_zres + vit_zres + gr_zres + ml_zres + at_zres + plot_layout(nrow = 4, ncol = 5, byrow = TRUE)
dev.off()

############## model predictions
# pulling predictions (may just want the values separately for now?)
# csv of predicted reservoirs and another for zoonotic?

## average predictions: Overall virus reservoir
virus_apreds=lapply(vbinary_brts,function(x) x$predict)
virus_apreds=do.call(rbind,virus_apreds)

## aggregate
virus_apreds=data.frame(aggregate(pred~species,data=virus_apreds,mean),
                        aggregate(cpred~species,data=virus_apreds,mean)['cpred'], ## holding wos constant
                        aggregate(dum_virus~species,data=virus_apreds,prod)["dum_virus"],
                        aggregate(dum_zvirus~species,data=virus_apreds,prod)["dum_zvirus"])

### type
# virus_apreds$type='PCR'

## average predictions: Zoonotic
zvirus_apreds=lapply(zbinary_brts,function(x) x$predict)
zvirus_apreds=do.call(rbind,zvirus_apreds)

## aggregate
zvirus_apreds=data.frame(aggregate(pred~species,data=zvirus_apreds,mean),
                         aggregate(cpred~species,data=zvirus_apreds,mean)['cpred'], ## holding wos constant
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


# pull relative importance
vinf <- lapply(vrichness_brts,function(x) x$rinf)

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

# calculate average rankings
ranks <- data_vinf %>%
  group_by(var) %>%
  summarize(avg = mean(rank),
            rse = std.error(rank),
            rvar = var(rank)) %>%
  ungroup() %>%
  arrange(avg)




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


# aggregate
agg = do.call(rbind, lapply(zbinary_brts,function(x) pdp_agg(x, "Synurbic")))

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
temp$predictor=temp["Synurbic"][,1]

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
