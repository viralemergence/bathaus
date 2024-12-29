# 00_descriptive stats
# Descriptive statistics and figure 1

# clear environment
rm(list=ls())
graphics.off()

# Load packages
library(tidyverse)
library(ggridges)

#### Descriptive stats 
# traits is the full dataset B4 cut off, has families in it
traits <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/synurbic and traits only.rds")

# add lowecase family names
traits$family <- tools::toTitleCase(tolower(traits$fam))

# density plots
traits %>% 
  filter(!fam %in% c("MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
  ggplot(aes(x = virus+1, y = fct_rev(family))) +
  geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
  #facet_wrap(~fam, ncol = 1) +
  scale_x_log10(limits = c(1, 1000)) +
  labs(x="Viral Richness (Log Scale)", y = " ", fill = "Roost Status") +
  #xlim(0, 1000) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.position = "none") +
  scale_fill_manual(labels = c("no","yes","missing"), values = c("#8470ff","#9DD866","#A0B1BA")) -> vrich

# zoonotic proportions
traits %>%
  mutate(zoo_prop = zvirus/virus) %>%
  mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop)) %>%
  filter(!fam %in% c("MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
  ggplot(aes(x = zoo_prop, y = fct_rev(family))) +
  geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
  #facet_wrap(~fam, ncol = 1) +
  #scale_x_log10() +
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x="Proportion Zoonotic Virus", y = " ", fill = "Roost Status") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text.y = element_blank()) +
  scale_fill_manual(labels = c("no","yes","missing"), values = c("#8470ff","#9DD866","#A0B1BA")) -> zprop

# multiplot
library(patchwork)
dens <- vrich + zprop

# save
png("/Volumes/BETKE 2021/bathaus/figs/figure 1e.png",width=5.5,height=6,units="in",res=600)
dens + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & 
  theme(legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),
        plot.tag = element_text(size = 8),
        legend.position = "bottom")
dev.off()

## Additional summary stats
# ranges of richness and zoonotic proportions
range(traits$virus)
mean(traits$virus)
median(traits$virus)

# number of known hosts
ct <- traits %>% select(species,virus,zvirus) %>%
  mutate(dum_virus = if_else(virus <= 0, 0, 1),
         dum_zvirus = if_else(zvirus <= 0, 0, 1),
         zoo_prop = zvirus/virus) 
sum(ct$dum_virus) # total - 341
sum(ct$dum_zvirus) # zoonotic - 205

# proportions 
sum(ct$dum_virus)/nrow(ct)
sum(ct$dum_zvirus)/nrow(ct)

# how many known hosts are anthropogenic roosting?
table(ct$dum_virus, ct$Synurbic, useNA = "ifany") # 227/341

# how many known zoonotic hosts are anthropogenic roosting
table(ct$dum_zvirus, ct$Synurbic, useNA = "ifany") # 175/266