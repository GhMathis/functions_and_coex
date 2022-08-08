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
library(glmmTMB)



localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "organized_data"

load(file.path(data_dir,"bioma.RData"))
load(file.path(data_dir,"coexistance.RData"))
load(file.path(data_dir,"cover.RData"))
load(file.path(data_dir,"par_and_height.RData"))
load(file.path(data_dir,"canopy_height_growth.RData"))
load(file.path(data_dir,"enzimatic_activity.RData"))
enz_act
########## OVERVIEW COEX DATA ##########

str(coexistance)
coexistance
ggplot(coexistance)+
  geom_point(aes(species, Omega, col = treatment))
ggplot(coexistance)+
  geom_point(aes(species, theta, col = treatment))
ggplot(coexistance)+
  geom_point(aes(Omega,theta , col = treatment))

main_data = coexistance
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
mod_sp_effect = lm(data = bioma_long, biomass ~ specie )
par(mfrow = c(2, 2))
plot(mod_sp_effect)
par(mfrow = c(1, 1))
summary(mod_sp_effect)
hist(mod_sp_effect$residuals)
### Fr and Hl happens to be in plot with a more biomass that other species
# => potential specie effect

head(main_data)
str(bioma)

main_data = full_join(main_data, bioma, by = c("treatment", "species"))

########## OVERVIEW COVER DATA ##########

str(cover)
CV <- function(cover){
  ###Community Variance per plot between the 3 sp
  var_per_plot =c()
  
  for (n_plot in 1:nrow(cover)){
    sp_vec = c(cover$species_1[n_plot], cover$species_2[n_plot], cover$species_3[n_plot])
    sp_vec = sp_vec[!is.na(sp_vec)]
    
    #variance of cover per plot
    var_per_plot = append(var_per_plot, sqrt(sd(as.matrix(cover[n_plot,sp_vec]))))
  }
  cover$var_per_plot = var_per_plot
  return(cover)
}



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

### DO the treatment interact have more impact on mono than on triplet on the
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

### Do species impact variance of the community cover ? and the treatment ?

#get the community cover variance 
cover<-CV(cover)

names(cover)

cover_long_var = cover[,-(6:20)] %>%
  pivot_longer(cols = species_1:species_3, names_to = NULL, values_to = "specie")
str(cover_long2)

mod_cover_var = lm(data = cover_long_var, var_per_plot ~ specie * treatment )

par(mfrow = c(2, 2))
plot(mod_cover_var)
par(mfrow = c(1, 1))
hist(residuals(mod_cover_var))

summary(mod_cover_var)

ggplot(cover_long_var)+
  geom_boxplot(aes(specie,var_per_plot, col=treatment))+
  ylab("Community cover variance")
### almost half of the species have an impact on the community cover variance
#(all of them increas the var between 1 and 2).
# Nitrogen also have a positive effect on it. But the interactions of nitrogen and 
# species have no impact.
str(main_data)
str(cover)
temp = cover[,c(1:2,28)]
temp$plot_ID = as.factor(temp$plot_ID)
main_data = full_join(main_data, temp, by = c("plot_ID", "treatment"))

########## OVERVIEW CANOPY HEIGHT DATA #############
str(canopy_height_growth)


canopy_growth_long = canopy_height_growth[,c(2,4:6)]%>%
  pivot_longer(cols= -(plot_ID:treatment),
               names_to = "coefs_names", values_to = "coefs")
  
ggplot(canopy_growth_long)+
  geom_boxplot(aes(treatment, coefs, col = coefs_names ))
mod_canopy_growth0=lm(data=canopy_height_growth, tang0_slope ~  treatment )
mod_canopy_growth38=lm(data=canopy_height_growth, tang38_slope ~  treatment )

par(mfrow = c(2, 2))
plot(mod_cover_var)
par(mfrow = c(1, 1))

summary(mod_canopy_growth0)
str(main_data)
str(canopy_height_growth)
temp = canopy_height_growth[,c(2,4:6)]
main_data = full_join(main_data, temp, by = c("plot_ID", "treatment"))
########### OVERVIEW ENZIMATICS ACTIVITIES DATA #############
beta_glu_df = subset(enz_act, enz_act$enz =="beta_glu")
phospha_df = subset(enz_act, enz_act$enz =="phospha")

str(beta_glu_df)
ggplot(beta_glu_df)+
  geom_boxplot(aes(cali_group, nytro_release))
ggplot(beta_glu_df)+
  geom_boxplot(aes(type, nytro_release))

ggplot(phospha_df)+
  geom_boxplot(aes(cali_group, nytro_release))
ggplot(phospha_df)+
  geom_boxplot(aes(type, nytro_release))

### calibration group have an strong impact that => random factor
# mono sems to realease less nytrophenol which seems logic.
temp_enz_df = data.frame(beta_glu = beta_glu_df$nytro_release,
                        phospha = phospha_df$nytro_release,
                        treatment = beta_glu_df$treatment,
                        plot_ID = beta_glu_df$plot_ID
)
ggplot(temp_enz_df)+
  geom_text(aes(beta_glu, phospha, col = treatment, label = plot_ID))+
  geom_abline(slope = 1, intercept = 0)+
  scale_y_continuous(limits = c(0,3))
#something went probably wrong for plot 5 (C and N), plot 3(N) and plot 4(N) for beta glu.
# MIGHT NEED TO BE REMOVE FROM THE DATA!
names(beta_glu_df)[1] = "nytro_release_B"
names(phospha_df)[1] = "nytro_release_P"
beta_glu_df[beta_glu_df$plot_ID %in% c(5,4,3) &
              beta_glu_df$treatment == "nitrogen", ]$nytro_release_B = NA
beta_glu_df[beta_glu_df$plot_ID == 5 &
              beta_glu_df$treatment == "control", ]$nytro_release_B = NA

str(main_data)
str(beta_glu_df)
temp = beta_glu_df[,c(1,3:4,8)]
main_data = full_join(main_data, temp, by = c("plot_ID", "treatment"))

temp = phospha_df[,c(1,3:4,8)]
main_data = full_join(main_data, temp, by = c("plot_ID", "treatment"))
names(main_data)[c(21,23)] =c("cali_grp_beta", "cali_grp_phospha")
head(main_data)
########### link functions and coexistance ##################
###### Normalisation of the data ######



norm.func <- function(df){
  df_norm = df
  for (col in 1:ncol(df)){
    max = max(df[,col], na.rm=T)
    min = min(df[,col], na.rm=T)
    df_norm[,col] = (df[,col] - min )/(max - min )
  }
  return(df_norm)
}
names(main_data)
main_data_norm = main_data
main_data_norm[,c(12,17:20,22,24:25)] = norm.func(main_data[,c(12,17:20,22,24:25)])

###### BIOMASSE ######

head(main_data_norm)
str(main_data_norm)
main_data_norm$feasibility = as.factor(main_data_norm$feasibility)
ggplot(main_data_norm)+
  geom_point(aes(Omega, biomass, col = treatment))+
  geom_smooth(aes(Omega, biomass), method = lm)+
  scale_x_log10()


ggplot(main_data_norm)+
  geom_point(aes(theta, biomass, col = treatment))
  geom_smooth(aes(theta, biomass), method = lm)+
  scale_x_log10()
ggplot(main_data_norm[main_data_norm$type!="mono",])+
  geom_boxplot(aes(feasibility, biomass, col = treatment))
hist(main_data_norm$Omega)

mod_biomass = glmmTMB(biomass ~(log(Omega) + log(theta) + differential )*treatment+
                      (1|species_1)+ (1|species_2) + (1|species_3),
                   data = main_data_norm)

par(mfrow = c(2, 2))
plot(mod_biomass)
par(mfrow = c(1, 1))
hist(residuals(mod_biomass))
summary(mod_biomass)
plot_model(mod_biomass, show.values = TRUE, value.offset = .3, 
           order.terms = c(1,6,2,7,3,8,4,9,5))

ggplot(main_data_norm,aes(differential, biomass, label= plot_ID))+
  geom_point()+
  geom_text()+
  geom_smooth(aes(differential, biomass), method = lm)+
  facet_wrap(~treatment)

##### GROWTH COEF #####

mod_tang0 = glmmTMB(tang0_slope ~(log(Omega) + log(theta) + differential )*treatment+
                      (1|species_1) + (1|species_2) + (1|species_3),
                   data = main_data_norm)

str(summary(mod_tang0))
df_coef=cbind(as.data.frame(summary(mod_tang0)$coefficients$cond[c(7,8),]),
              func = "tang0")
df_coef = df_coef[1,]
plot_model(mod_tang0, show.values = TRUE, value.offset = .3, 
           order.terms = c(1,6,2,7,3,8,4,9,5))

mod_max_heigth = glmmTMB(max_canopy_heigth ~(log(Omega) + log(theta) + differential )*treatment+
                           (1|species_1) + (1|species_2) + (1|species_3),
                 data = main_data_norm)

summary(mod_max_heigth)
df_coef=rbind(df_coef,cbind(as.data.frame(summary(mod_max_heigth)$coefficients$cond[c(2,4),]),
      func = "max_heigth"))
plot_model(mod_max_heigth, show.values = TRUE, value.offset = .3,
           order.terms = c(1,6,2,7,3,8,4,9,5))
head(main_data_norm)
mod_max_time = glmmTMB(max_canopy_time ~(log(Omega) + log(theta) + differential )*treatment+
                           (1|species_1) + (1|species_2) + (1|species_3),
                         data = main_data_norm)

summary(mod_max_time)
df_coef=rbind(df_coef,cbind(as.data.frame(summary(mod_max_time)$coefficients$cond[c(2,4),]),
                            func = "max_time"))

plot_model(mod_max_heigth, show.values = TRUE, value.offset = .3,
           order.terms = c(1,6,2,7,3,8,4,9,5))


ggplot(main_data_norm)+
  geom_point(aes(Omega,var_per_plot))+
  facet_wrap(~treatment)
mod_var_cover = glmmTMB(var_per_plot ~(log(Omega) + log(theta) + differential )*treatment+
                          (1|species_1) + (1|species_2) + (1|species_3),
                 data = main_data_norm)

summary(mod_var_cover)
plot_model(mod_var_cover, show.values = TRUE, value.offset = .3,
           order.terms = c(1,6,2,7,3,8,4,9,5))



mod_beta = glmmTMB(nytro_release_B ~(log(Omega) + log(theta) + differential )*treatment+
                       (1|species_1) + (1|species_2) + (1|species_3)+(1|cali_grp_beta),
                     data = main_data_norm)

summary(mod_beta)
df_coef=rbind(df_coef,cbind(as.data.frame(summary(mod_beta)$coefficients$cond[4:5,]),
                            func = "beta"))
df_coef = df_coef[1:6,]
plot_model(mod_beta, show.values = TRUE, value.offset = .3,
           order.terms = c(1,6,2,7,3,8,4,9,5))
mod_phospha = glmmTMB(nytro_release_P ~ (log(Omega) + log(theta) + differential ) * treatment+
                        (1|cali_grp_phospha),
                        data = main_data_norm)
ggplot(main_data_norm)+
  geom_point(aes(differential, nytro_release_P))+
  facet_wrap(~treatment)
  scale_x_log10()
summary(mod_phospha)
plot_model(mod_phospha, show.values = TRUE, value.offset = .3,
           order.terms = c(1,6,2,7,3,8,4,9,5))
