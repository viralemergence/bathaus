# 00_descriptive stats
# Descriptive statistics and figure 1

rm(list=ls())
graphics.off()

library(tidyverse)
library(ggridges)

#### Descriptive stats
# traits is the full dataset B4 cut off, has families in it
traits <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/synurbic and traits only.rds")

# add lowecase family names
traits$family <- tools::toTitleCase(tolower(traits$fam))

# density plots
traits %>% filter(!fam %in% c("MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
  ggplot(aes(x = virus+1, y = fct_rev(Family))) +
  geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
  #facet_wrap(~fam, ncol = 1) +
  scale_x_log10(limits = c(1, 1000)) +
  labs(x="Viral Richness (Log Scale)", y = " ", fill = "Roost Status") +
  #xlim(0, 1000) +
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
  filter(!fam %in% c("MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
  ggplot(aes(x = zoo_prop, y = fct_rev(Family))) +
  geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
  #facet_wrap(~fam, ncol = 1) +
  #scale_x_log10() +
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
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

# try with zeros turned to NA for frequencies
na_points <- points %>% mutate(clean = na_if(Frequency, 0))

# I could pull out to match the other two plots
#points[!points$Var1 %in% c('NYCTERIDAE',"MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE"), ]

ggplot(na_points, aes(x = `Roost Status`, y = fct_rev(Family))) +
  geom_point(aes(size = clean, color = `Roost Status`)) +
  labs(x="Anthropogenic \nRoosting", y = " ", color = "Anthropogenic \nRoosting") +
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
  scale_size(range = c(0.5, 6), name = "Frequency") +
  guides(color = guide_legend(override.aes = list(size=3))) -> balloon1

library(patchwork)
dens <- vrich + zprop

# put together
fig1 <- (balloon1 + plot_layout(widths = c(15, 85), guide = "collect") & 
           theme(legend.position = "right",
                 #legend.box = "vertical",
                 legend.text = element_text(size = 6), 
                 legend.title = element_text(size = 6),
                 legend.margin=margin(l=-7))) + 
  (vrich + zprop)

png("~/Desktop/Bats and Viruses/bathaus/figs/figure 1.png",width=6,height=3.5,units="in",res=600)
fig1 + plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 7))
dev.off()

# # lighthouse slides panels B1 and B2 only
# # density plots
# traits %>% filter(!fam %in% c('NYCTERIDAE',"MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
#   ggplot(aes(x = virus+1, y = family)) +
#   geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
#   #facet_wrap(~fam, ncol = 1) +
#   scale_x_log10() +
#   labs(x="Log Viral Richness", y = " ", fill = "Anthropogenic Roost") +
#   theme_bw() +
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         axis.text = element_text(size = 6),
#         axis.title = element_text(size = 7),
#         legend.position = "none") +
#   scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) -> vrich2

# # zoonotic proportions
# traits %>%
#   mutate(zoo_prop = zvirus/virus) %>%
#   mutate(zoo_prop = ifelse(is.nan(zoo_prop), 0, zoo_prop)) %>%
#   filter(!fam %in% c('NYCTERIDAE',"MYSTACINIDAE","CRASEONYCTERIDAE","FURIPTERIDAE","MYZOPODIDAE", "NOCTILIONIDAE")) %>%
#   ggplot(aes(x = zoo_prop, y = family)) +
#   geom_density_ridges(aes(fill = Synurbic), alpha = 0.6) +
#   #facet_wrap(~fam, ncol = 1) +
#   #scale_x_log10() +
#   scale_x_continuous(labels = scales::percent) +
#   labs(x="Proportion Zoonotic Virus", y = " ", fill = "Anthropogenic Roost") +
#   theme_bw() +
#   theme(panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         axis.text = element_text(size = 6),
#         axis.title = element_text(size = 7),
#         legend.position = "right",
#         axis.text.y = element_blank()) +
#   scale_fill_manual(labels = c("no","yes","NA"), values = c("#8470ff","#9DD866","#A0B1BA")) -> zprop2
# 
# png("/Volumes/BETKE 2021/bathaus/figs/Lighthouse fig 1.png",width=5.5,height=3.5,units="in",res=600)
# (vrich2 + zprop2) + theme(legend.position = "right",
#                         #legend.box = "vertical",
#                         legend.text = element_text(size = 5), 
#                         legend.title = element_text(size = 6),
#                         legend.key.size = unit(0.3, "cm"))
# dev.off()

## Additional summary stats
data <- readRDS("~/Desktop/Bats and Viruses/bathaus/flat files/cleaned dataset 30 cutoff.rds")

## ranges of richness and zoonotic proportions
range(data$virus)
mean(data$virus)
median(data$virus)

# number of known hosts
data <- data %>%
  mutate(dum_virus = if_else(virus <= 0, 0, 1),
         dum_zvirus = if_else(zvirus <= 0, 0, 1)) 
sum(data$dum_virus) # total - 381
sum(data$dum_zvirus) # zoonotic - 266

# proportions 
sum(data$dum_virus)/nrow(data)
sum(data$dum_zvirus)/nrow(data)

# how many known hosts are anthropogenic roosting?
table(data$dum_virus, data$Synurbic, useNA = "ifany") # 243/379

# how many known zoonotic hosts are anthropogenic roosting
table(data$dum_zvirus, data$Synurbic, useNA = "ifany") # 175/266

# quick boxplots?
ggplot(data) +
  geom_boxplot(aes(x=factor(Synurbic),y=virus+1, fill=factor(Synurbic))) +
  scale_y_log10() +
  theme_bw()

# filtered 
data %>%
  filter(virus >= 1) %>% # knowns 
  ggplot() +
  geom_boxplot(aes(x=factor(Synurbic),y=virus,fill=factor(Synurbic))) +
  scale_y_log10() +
  theme_bw()

# zoonotic virus
ggplot(data) +
  geom_boxplot(aes(x=factor(Synurbic),y=log(zvirus+1))) +
  theme_bw()

# filtered to knowns
data %>%
  filter(zvirus >= 1) %>% # knowns 
  ggplot() +
  geom_boxplot(aes(x=factor(Synurbic),y=zvirus,fill=factor(Synurbic))) +
  scale_y_log10() +
  theme_bw()

## Can you break down the values of the plot by family - might want to get 
# so this would possibly be average richness for each synurbic type - not very descriptive
traits %>% 
  filter(factor(Synurbic) == 1) %>%
  group_by(fam) %>%
  summarize(avg = mean(virus),
            med = median(virus)) %>%
  ungroup()

traits %>% 
  filter(factor(Synurbic) == 0) %>%
  group_by(fam) %>%
  summarize(avg = mean(virus),
            med = median(virus)) %>%
  ungroup()