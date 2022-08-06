library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)
library(lmerTest)
library(ade4)
library(partR2)
library(tidyverse)
library(sjPlot)



localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "organized_data"

load(file.path(data_dir,"bioma.RData"))
load(file.path(data_dir,"coexistance.RData"))
load(file.path(data_dir,"cover.RData"))



########## OVERVIEW COEX DATA ##########

str(coexistance)
coexistance
ggplot(coexistance)+
  geom_point(aes(species, Omega, col = treatment))
ggplot(coexistance)+
  geom_point(aes(species, theta, col = treatment))
ggplot(coexistance)+
  geom_point(aes(Omega,theta , col = treatment))

########## OVERVIEW BIOMASSE DATA ##########
head(bioma)
str(bioma)

bioma_long = bioma %>%
  pivot_longer(cols = species_1:species_3, names_to = "sp_number", values_to = "specie" )
head(bioma_long)
ggplot(bioma_long)+
  geom_point(aes(plot_ID, biomass, col = treatment))+
  facet_wrap(~type)
ggplot(subset(bioma_long, bioma_long$type == "triplet"))+
  geom_boxplot(aes(specie, biomass))+
  facet_wrap(~treatment)
ggplot(bioma_long)+
  geom_point(aes(plot_ID, biomass, col = treatment))+
  facet_wrap(~type)

str(bioma_long)
mod_sp_effect = lm(data = bioma_long, biomass ~ specie)
par(mfrow = c(2, 2))
plot(mod_sp_effect)
par(mfrow = c(1, 1))
summary(mod_sp_effect)
hist(mod_sp_effect$residuals)
### Fr and Hl happens to be in plot with a more biomass that other species
# => potential specie effect


######## biomass and coex ##########

str(coexistance)
str(bioma)

data_biomass = full_join(coexistance, bioma, by = c("treatment", "species"))
str(data_biomass)
ggplot(data_biomass)+
  geom_point(aes(Omega, biomass, col = treatment))+
  scale_x_log10()
ggplot(data_biomass)+
  geom_point(aes(theta, biomass, col = treatment))
  

########## OVERVIEW COVER DATA ##########

str(cover)
CV <- function(cover){
  ###Community Variance variance per plot between the 3 sp
  var_per_plot =c()
  
  for (n_plot in 1:nrow(cover)){
    sp_vec = c(cover$species_1[n_plot], cover$species_2[n_plot], cover$species_3[n_plot])
    sp_vec = sp_vec[!is.na(sp_vec)]
    var_per_plot = append(var_per_plot, sqrt(sd(as.matrix(cover[n_plot,sp_vec]))))
  }
  cover$var_per_plot = var_per_plot
  return(cover)
}

cover<-CV(cover)


cover_long = cover %>%
  pivot_longer(cols = Be:legumes, names_to = "specie_for_cover", values_to = "perc_cover")


ggplot(cover_long)+
  geom_boxplot(aes(specie_for_cover, perc_cover))
str(cover_long)
mod_cover = lm(data = cover_long[!(cover_long$specie_for_cover  %in%
                                   c("grasses", "herbs","legumes")),]
               ,perc_cover ~ specie_for_cover)

par(mfrow = c(2, 2))
plot(mod_cover)
par(mfrow = c(1, 1))
hist(mod_cover$residuals)
summary(mod_cover)
coef(mod_cover)
plot_model(mod_cover,show.values = TRUE, value.offset = .3, sort.est = TRUE)

### DO the treatment interact have more impact on mon than on triplet on the
# the percentage cover ?
mod_cover2 = lm(data = cover_long[!(cover_long$specie_for_cover  %in%
                                     c("grasses", "herbs","legumes")),]
               ,perc_cover ~ treatment * type )

par(mfrow = c(2, 2))
plot(mod_cover2)
par(mfrow = c(1, 1))
hist(mod_cover2$residuals)
summary(mod_cover2)
coef(mod_cover2)
str(cover_long)

ggplot(cover_long[!(cover_long$specie_for_cover  %in%
                      c("grasses", "herbs","legumes")),])+
  geom_boxplot(aes(treatment,perc_cover, col=type))+
  geom_point(aes(treatment,perc_cover, col=type))



