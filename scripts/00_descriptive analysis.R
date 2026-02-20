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
traits <- readRDS("~/Desktop/bathaus/flat files/synurbic and traits only.rds")

# add lowecase family names
traits$family <- tools::toTitleCase(tolower(traits$fam))

# density plots
traits %>% 
  filter(!fam %in% c("MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
  ggplot(aes(x = vfam+1, y = fct_rev(family))) +
  geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
  #facet_wrap(~fam, ncol = 1) +
  scale_x_log10(limits = c(1, 150)) +
  labs(x="Virus Family Richness (Log Scale)", y = " ", fill = "Roost Status") +
  #xlim(0, 1000) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        legend.position = "none") +
  scale_fill_manual(labels = c("no","yes","missing"), values = c("#8470ff","#9DD866","#A0B1BA")) -> vrich

# zoonotic family richness
traits %>%
  # mutate(zoo_prop = zvirus/virus) %>%
  # mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop)) %>%
  filter(!fam %in% c("MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
  ggplot(aes(x = zfam+1, y = fct_rev(family))) +
  geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
  #facet_wrap(~fam, ncol = 1) +
  #scale_x_log10() +
  #scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_x_log10(limits = c(1, 150)) +
  labs(x="Zoonotic Virus Famiy Richness (Log Scale)", y = " ", fill = "Roost Status") +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        axis.text.y = element_blank()) +
  scale_fill_manual(labels = c("no","yes","missing"), values = c("#8470ff","#9DD866","#A0B1BA")) -> zrich

# multiplot
library(patchwork)
dens <- vrich + zrich

# save
png("/Users/brianabetke/Desktop/bathaus/figs/figure 1.png",width=5.5,height=6,units="in",res=600)
dens + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A') & 
  theme(legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),
        plot.tag = element_text(size = 8),
        legend.position = "bottom")
dev.off()

## Additional summary stats
# ranges of richness and zoonotic proportions
range(traits$vfam)
mean(traits$vfam)
median(traits$vfam)

# number of known hosts
# ct <- traits %>% select(species,virus,zvirus) %>%
#   mutate(dum_virus = if_else(virus <= 0, 0, 1),
#          dum_zvirus = if_else(zvirus <= 0, 0, 1),
#          zoo_prop = zvirus/virus) 
sum(traits$dum_virus) # total - 367
sum(traits$dum_zvirus) # zoonotic - 222

# proportions 
sum(traits$dum_virus)/nrow(traits) # 0.2869429
sum(traits$dum_zvirus)/nrow(traits) #  0.1735731

# how many known hosts are anthropogenic roosting?
table(traits$dum_virus, traits$Synurbic, useNA = "ifany") # 243/367

# how many known zoonotic hosts are anthropogenic roosting
table(traits$dum_zvirus, traits$Synurbic, useNA = "ifany") # 151/222
