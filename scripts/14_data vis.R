# data vis - descriptive and brts?

#clean envrionment
rm(list=ls())
graphics.off()

# test code for pulling out all the relative importance info from every model
library(tidyverse)
library(plotrix)
library(patchwork)
library(gbm)
#library(ggrepel)

#### Descriptive stats
library(ggridges)
data <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/cleaned dataset 30 cutoff.rds")
# traits is the full dataset B4 cut off, has families in it
traits <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/synurbic and traits only.rds")

# add lowecase family names
traits$family <- tools::toTitleCase(tolower(traits$fam))

# density plots
traits %>% filter(!fam %in% c('NYCTERIDAE',"MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
  ggplot(aes(x = virus+1, y = family)) +
  geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
  #facet_wrap(~fam, ncol = 1) +
  scale_x_log10() +
  labs(x="Log Viral Richness", y = " ", fill = "Roost Status") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.position = "none") +
scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) -> vrich

# zoonotic proportions
traits %>%
  mutate(zoo_prop = zvirus/virus) %>%
  mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop)) %>%
  filter(!fam %in% c('NYCTERIDAE',"MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
  ggplot(aes(x = zoo_prop, y = family)) +
  geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
  #facet_wrap(~fam, ncol = 1) +
  #scale_x_log10() +
  scale_x_continuous(labels = scales::percent) +
  labs(x="Proportion Zoonotic Virus", y = " ", fill = "Roost Status") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7),
        legend.position = "none",
        axis.text.y = element_blank()) +
  scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) -> zprop

#balloon plot
points <- data.frame(table(traits$family,traits$Synurbic, useNA = "ifany"))
colnames(points) <- c("Family", "Roost Status", "Frequency")
  
# I could pull out to match the other two plots
#points[!points$Var1 %in% c('NYCTERIDAE',"MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE"), ]

ggplot(points, aes(x = `Roost Status`, y = Family)) +
  geom_point(aes(size = Frequency, color = `Roost Status`)) +
  labs(x="Roost Status", y = " ", fill = "Roost Status") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7)) +
        # legend.position = "top",
        # legend.box = "vertical",
        # legend.margin=margin()) +
  scale_color_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_x_discrete(labels = c('no',"yes","NA")) +
  scale_size(range = c(0.5, 6)) +
  guides(color = guide_legend(override.aes = list(size=3))) -> balloon1

# # I could pull out to match the other two plots
# trim <- points[!points$Var1 %in% c('NYCTERIDAE',"MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE","NOCTILIONIDAE"), ]
# 
# ggplot(trim, aes(x = Var2, y = Var1)) +
#   geom_point(aes(size = Freq, color = Var2)) +
#   labs(x="Roost Status", y = " ", fill = "Roost Status") +
#   theme_bw() +
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         legend.position = "top", 
#         legend.box="vertical",
#         legend.margin=margin()) +
#   scale_color_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
#   scale_size(range = c(1, 8))-> balloon2 # to match other panels, turns out this doesn't work the way I thought....

ggplot(points, aes(x = `Roost Status`, y = Family)) +
  geom_point(aes(size = Frequency, color = `Roost Status`)) +
  labs(x=NULL, y =NULL, fill = "Roost Status") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 7)) +
  # legend.position = "top",
  # legend.box = "vertical",
  # legend.margin=margin()) +
  scale_color_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) +
  scale_x_discrete(labels = c('no',"yes","NA")) +
  scale_size(range = c(1, 5)) +
  guides(color = guide_legend(override.aes = list(size=3))) -> balloon

dens <- vrich + zprop

# put together
fig1 <- (balloon1 + plot_layout(widths = c(15, 85), guide = "collect") & 
    theme(legend.position = "right",
          #legend.box = "vertical",
          legend.text = element_text(size = 6), 
          legend.title = element_text(size = 6),
          legend.margin=margin(l=-7))) + 
  ((vrich + zprop) + plot_layout(tag_level = "new"))

(balloon1 + dens) + plot_layout(widths = c(15, 85), guide = "collect") & 
  theme(legend.position = "bottom",
        #legend.box = "vertical",
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),
        legend.margin=margin(b = 0, unit='cm'))

alt <- (balloon + 
    plot_layout(heights = c(15, 85), guide = "collect") & 
    theme(legend.position = "bottom",
          legend.text = element_text(size = 6), 
          legend.title = element_text(size = 7),
          legend.margin=margin(l=-7))) + 
  ((vrich + zprop) + plot_layout(tag_level = "new"))

(balloon1 + vrich) + plot_layout(widths = c(15, 85), guide = "collect") & 
    theme(legend.position = "right",
          #legend.box = "vertical",
          legend.text = element_text(size = 6), 
          legend.title = element_text(size = 6),
          legend.margin=margin(l=-7)) + 
   (zprop + plot_layout(tag_level = "new"))

png("~/Desktop/Bats and Viruses/bathaus/figs/Figure 1.png",width=6.5,height=3.5,units="in",res=600)
fig1 + plot_annotation(tag_levels = c('A','1')) & theme(plot.tag = element_text(size = 7))
dev.off()

png("~/Desktop/Bats and Viruses/bathaus/figs/Figure 1 Alt.png",width=6,height=6,units="in",res=600)
alt + plot_annotation(tag_levels = list(c('A','B', NULL))) & theme(plot.tag = element_text(size = 7), plot.tag.position = "topleft")
dev.off()

###################### BRT Results
# read in the datasets 
fvirus_brts <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/fvirus brts.rds")
fzvirus_brts <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/fzvirus brts.rds")

# full brts
virus_brts <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/virus brts.rds")
zvirus_brts <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/zvirus brts.rds")

# Read in rds files for model outputs
vrichness_brts <- readRDS("/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/virus with brts.rds")
no_vrichness_brts <- readRDS("/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/virus without brts.rds")
zoo_prop_brts <- readRDS("/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/zoo_prop with brts.rds")
no_zoo_prop_brts <- readRDS("/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/zoo_prop without brts.rds")
vbinary_brts <- readRDS("/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/dum_virus with brts.rds")
no_vbinary_brts <- readRDS("/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/dum_virus without brts.rds")
zbinary_brts <- readRDS("/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/dum_zvirus with brts.rds")
no_zbinary_brts <- readRDS("/Users/brianabetke/Desktop/Bats and Viruses/bathaus/flat files/dum_zvirus without brts.rds")

################### Variable Importance Plots
# Pull all the relative importance into a dataframe, get the mean, sd, and variation.
# Then create a plot similar to the one I made for the variants 
vinfPlot <- function(data_name, df_name, fig_name, bar_color){
  
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
  
  # gather Test AUC
  tauc <- sapply(data_name, function(x) x$testAUC)
  mean_tauc <- mean(tauc)
  se_tauc <- std.error(tauc)
  test_auc <- data.frame(avg = mean_tauc, 
                         se = se_tauc)
  # vinf=lapply(data_name,function(x) x$rinf)
  # data_vinf=do.call(rbind,vinf)
  # ## from brtvis aggregate mean, SE, var
  # # creates a dataset from list with 10 unique splits
  # vdata=data.frame(aggregate(rel.inf~var,data=data_vinf,mean),
  #                       aggregate(rel.inf~var,data=data_vinf,st.error)["rel.inf"],
  #                       aggregate(rel.inf~var,data=data_vinf,var)["rel.inf"])
  # names(vdata)=c("var","rel.inf","rse","rvar")
  # df_name=vdata[order(vdata$rel.inf,decreasing=T),]
  
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
  
  fig_name <- ggplot(df_name, aes(x = reorder(var, -avg), y = avg, color = ifelse(var == "Anthropogenic Roost", "Anthro", "Shade"))) + 
    #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
    #geom_bar(stat = "identity") +
    #geom_errorbar(aes(ymin = avg-rse, ymax = avg+rse)) +
    geom_pointrange(aes(ymin = avg-rse, ymax = avg+rse)) +
    #coord_flip() +
    theme_bw() +
    theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
          legend.position = "none") +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(size = 10)) + 
    scale_color_manual(breaks = c("Anthro","Shade"), values = c("red", "purple"))
  
  # return a list with that dataset of rel.inf and figure
  return(list(df_name, test_auc, fig_name))
}

comp_fig <- vinfPlot(vrichness_brts, compdf, compfig, "grey")

# run function with dataset
# comp test
comp_fig <- vinfPlot(comp_brts, compdf, compfig, "grey")
c <- comp_fig[[3]] + 
  labs(x = " ", y = "Relative Importance") +
  theme(axis.title.y = element_text(size = 14, hjust = -1.5))

# pcr test
pcr_fig <- vinfPlot(pcr_brts, pcrdf, pcrfig, "palegreen3")
p <- pcr_fig[[3]] + 
  labs(x = "Trait Variable", y = " ") +
  theme(axis.title.x = element_text(size = 14))

# put the figures together
png("variableinf.png",width=15,height=10,units="in",res=600)
c / p
dev.off()

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
h_patch <- fvirus / fzvirus & ylab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

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

################# may want to calculate feature rankings first - get the top most important predictors in the model







################# going to want code here that makes pdps for all models
# work on adapting the synurbat code to handle aggregated data - will aggregating result in a file that looks the same but just of averages?
# Need to look at hantavirus code 
# the thing is calculating the marginal effects across models do I actually show averaged y axis?




## Partial dependence plots for continuous variables 
make_pdp_cont <- function(model, predictor, var_name, pcolor = FALSE) {
  
  # return grid
  vals <- plot.gbm(model[["mod"]], i.var = predictor, type = "response", return.grid = TRUE)
  
  # data for hist
  yrange = range(vals$y, na.rm = TRUE)
  
  # pull histogram values
  hi=hist(model[["testdata"]][[predictor]],breaks=30,plot=F)
  hi=with(hi,data.frame(breaks[1:(length(breaks)-1)],counts))
  names(hi)=c("mids","counts")
  
  if(pcolor == FALSE){
    
    # plot
    ggplot() + 
      geom_line(data = vals, aes(x = !!sym(predictor), y = y)) +
      geom_segment(data=hi,inherit.aes=F,
                   aes(x=mids,xend=mids,
                       y=yrange[1],yend=plotrix::rescale(counts,yrange)),
                   size=2,colour="darkgrey",alpha=0.50) +
      labs(x = var_name, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    
  }else{
    
    ggplot() + 
      geom_line(data = vals, aes(x = !!sym(predictor), y = y)) +
      geom_segment(data=hi,inherit.aes=F,
                   aes(x=mids,xend=mids,
                       y=yrange[1],yend=plotrix::rescale(counts,yrange)),
                   size=2,colour="Orange",alpha=0.40) +
      labs(x = var_name, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
  }
  
}

## Function for factor pdp plots
make_pdp_fact <- function(model, predictor, var_name, pcolor = FALSE) {
  
  # return grid
  vals <- plot.gbm(model[["mod"]], i.var = predictor, type = "response", return.grid = TRUE)
  
  # data for hist
  yrange = range(vals$y, na.rm = TRUE)
  
  # pull counts for color
  df_cat <- as.data.frame(table(noNA[["testdata"]][[predictor]]))
  
  # fix y axis point
  df_cat$ymin <- yrange[1]-0.01
  
  if(pcolor == FALSE){ #greys for initial model
    
    ggplot() +
      geom_point(data = vals, size= 2, shape=15, aes(category, y)) +
      geom_point(data = df_cat, aes(Var1, ymin, color = Freq)) +
      scale_color_continuous(high = "#636363", low = "#D9D9D9", guide = "none") +
      labs(x = var_name, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    
  }else{ #different color for pseudoab variables
    
    ggplot() +
      geom_point(data = vals, size= 2, shape=15, aes(category, y)) +
      geom_point(data = df_cat, aes(Var1, ymin, color = Freq)) +
      scale_color_continuous(high = "#8C2D04", low = "#FEE6CE", guide = "none") +
      labs(x = var_name, y = "Marginal Effect") +
      theme_bw() +
      theme(axis.text=element_text(size=6),
            axis.title=element_text(size=7)) +
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank())
    
  }
  
}

## Plot partial dependence
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


############### not sure if this is the order I want to put this here but take performance values and calculate sig. changes in performance metrics
# Cohen's D





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

# TESTING CODE
fig_name <- ggplot(g, aes(x = reorder(var, avg), y = avg, fill = ifelse(var == "Synurbic", "Anthro", "Shade"))) + 
  #geom_crossbar(aes(ymin = avg-rse, ymax = avg+rse), alpha = 0.5) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = avg-rse, ymax = avg+rse)) +
  #geom_pointrange(aes(ymin = avg-rse, ymax = avg+rse)) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6, hjust = 1),
        legend.position = "none") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(size = 10)) + 
  scale_fill_manual(breaks = c("Anthro","Shade"), values = c("red", "purple")) +
  labs(x = NULL, y = "sqrt avg variable influence")



fig_name + scale_y_sqrt()

