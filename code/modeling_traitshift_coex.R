library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)
library(lmerTest)
library(ade4)
library(partR2)
library(tidyverse)

localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)

data_dir = "organized_data"
code_dir ="code"

##### load the data
load(file.path(data_dir,"traits_shift_data.RData"))

###################################################
    ########## BIOMASSE ANNALYSYS ##########
###################################################

head(traits_shift_data)

##### models CWM #####
biomass_data = traits_shift_data[traits_shift_data$trait == "biomass",]
hist(biomass_data$CWM,breaks =48 )
hist(biomass_data$delta_ab_shift,breaks =48 )
hist(biomass_data$delta_intra_shift,breaks =48 )
mod_CWM_bioma <- lm(formula = CWM ~ Omega + 
                feasibility + 
                indirect_effects +
                log(year_extinction), data = biomass_data)
summary(mod_CWM_bioma)
coef_CWM_bioma = summary(mod_CWM_bioma)$coefficients

hist(summary(mod_CWM_bioma)$residuals, breaks = 48)
hist(biomass_data$CWM)
ggplot(biomass_data)+
  geom_point(aes(Omega, CWM))+
  geom_abline(intercept = coef_CWM_bioma[1,1], slope = coef_CWM_bioma[2,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_boxplot(aes(as.factor(feasibility), CWM))+
  geom_abline(intercept = coef_CWM_bioma[1,1], slope = coef_CWM_bioma[3,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_point(aes(indirect_effects, CWM))+
  geom_abline(intercept = coef_CWM_bioma[1,1], slope = coef_CWM_bioma[4,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_point(aes(year_extinction , CWM))+
  geom_abline(intercept = coef_CWM_bioma[1,1], slope = coef_CWM_bioma[5,1],
              linetype =2 , color = "blue")+
  scale_x_log10()


###################################################
  ########## beta glucosidase ANNALYSYS ##########
###################################################

head(traits_shift_data)

##### models CWM #####
beta_glu_data = traits_shift_data[traits_shift_data$trait == "beta_glu",]
hist(beta_glu_data$CWM,breaks =48 )
hist(beta_glu_data$delta_ab_shift,breaks =48 )
hist(beta_glu_data$delta_intra_shift,breaks =48 )
mod_CWM_betaG <- lm(formula = CWM ~ Omega + 
                feasibility + 
                indirect_effects +
                log(year_extinction), data = beta_glu_data)
summary(mod_CWM_betaG)
coef_CWM_betaG = summary(mod_CWM_betaG)$coefficients

hist(summary(mod_CWM_betaG)$residuals, breaks = 48)
hist(beta_glu_data$CWM)
ggplot(beta_glu_data)+
  geom_point(aes(Omega, CWM))+
  geom_abline(intercept = coef_CWM_betaG[1,1], slope = coef_CWM_betaG[2,1],
              linetype =2 , color = "blue")

ggplot(beta_glu_data)+
  geom_boxplot(aes(as.factor(feasibility), CWM))+
  geom_abline(intercept = coef_CWM_betaG[1,1], slope = coef_CWM_betaG[3,1],
              linetype =2 , color = "blue")

ggplot(beta_glu_data)+
  geom_point(aes(indirect_effects, CWM))+
  geom_abline(intercept = coef_CWM_betaG[1,1], slope = coef_CWM_betaG[4,1],
              linetype =2 , color = "blue")

ggplot(beta_glu_data)+
  geom_point(aes(year_extinction , CWM))+
  geom_abline(intercept = coef_CWM_betaG[1,1], slope = coef_CWM_betaG[5,1],
              linetype =2 , color = "blue")+
  scale_x_log10()

###################################################
########## phosphatase ANNALYSYS ##########
###################################################

head(traits_shift_data)
traits_shift_data$trait
##### models CWM #####
phospha_data = traits_shift_data[traits_shift_data$trait == "beta_glu",]
hist(phospha_data$CWM,breaks =48 )
hist(phospha_data$delta_ab_shift,breaks =48 )
hist(phospha_data$delta_intra_shift,breaks =48 )
mod_CWM_phospha <- lm(formula = CWM ~ Omega + 
                feasibility + 
                indirect_effects +
                log(year_extinction), data = phospha_data)
summary(mod_CWM_phospha)
coef_CWM_phospha = summary(mod_CWM_phospha)$coefficients

hist(summary(mod_CWM_phospha)$residuals, breaks = 48)
hist(phospha_data$CWM)
ggplot(phospha_data)+
  geom_point(aes(Omega, CWM))+
  geom_abline(intercept = coef_CWM_phospha[1,1], slope = coef_CWM_phospha[2,1],
              linetype =2 , color = "blue")

ggplot(phospha_data)+
  geom_boxplot(aes(as.factor(feasibility), CWM))+
  geom_abline(intercept = coef_CWM_phospha[1,1], slope = coef_CWM_phospha[3,1],
              linetype =2 , color = "blue")

ggplot(phospha_data)+
  geom_point(aes(indirect_effects, CWM))+
  geom_abline(intercept = coef_CWM_phospha[1,1], slope = coef_CWM_phospha[4,1],
              linetype =2 , color = "blue")

ggplot(phospha_data)+
  geom_point(aes(year_extinction , CWM))+
  geom_abline(intercept = coef_CWM_phospha[1,1], slope = coef_CWM_phospha[5,1],
              linetype =2 , color = "blue")+
  scale_x_log10()

