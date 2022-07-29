library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(emmeans)
library(lmerTest)
library(ade4)
library(partR2)
library(tidyverse)
###################################################
    ########## BIOMASSE ANNALYSYS ##########
###################################################

localDir="~/Fac/Cesure2/Plant species coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "organized_data"

load(file.path(data_dir,"bioma.RData"))

load(file.path(data_dir,"interactions_triplet.RData"))
load(file.path(data_dir,"cover.RData"))

##### shape data #####
names(interactions)[1] = "nitorgen" 
interactions$nitorgen = as.factor(rep(0,nrow(interactions )))

cover$species = apply( cover[,3:5], 1, 
                      function(x) paste(x[!is.na(x)], collapse = "_"))
cover$type = ifelse(as.numeric(as.character(cover$plot)) <= 48 ,
                   "triplet", "mono")
cover_long = cover %>%
  pivot_longer(cols = Be:legumes, names_to = "species_cover", values_to = "perc_cover")

str(cover_long)

str(cover)
hist(cover_long$perc_cover)
ggplot(cover_long)+
  geom_point(aes(species_cover,perc_cover))

head(interactions)
bioma$species = apply( bioma[,4:6], 1, 
                           function(x) paste(x[!is.na(x)], collapse = "_"))
str(bioma)
##### functions #####
sp_names = names(cover[1,6:17])
strategies <- setNames(c("slow", "fast", "slow", "slow", "fast", "fast", "slow",
                         "slow", "fast", "fast", "slow", "fast"),sp_names)

biomass.per.sp <- function(bioma, cover){
  ### function to find the biomass per species with the biomass data
  ### and the cover data.
  sp_cover = sum(cover[,6:17],na.rm = T)
  total_cover = cover$total.cover
  alien_cover = total_cover - sp_cover # cover of sp that aren't in the experiment (if needed)
  
  # cover weighted biomass calcul 
  bioma_weighted = (cover[,6:17]/total_cover) * bioma$biomasse
  
  return(as.matrix(bioma_weighted))
}
str(bioma_long)

biomass.shift <- function(bioma_long){
  ### fonction that find the biomass shift for each species
  ### bioma_long is the biomass data frame in a long format
  biomass_shift = c()
  for  (n in 1 : nrow(bioma_long)){
    sp = bioma_long[n, ]$species
    nitro = bioma_long[n, ]$nitrogen
    
    biomass_mono = subset(bioma_long, bioma_long$type == "mono" & 
                           bioma_long$species_1 == sp & 
                           bioma_long$species == sp &
                           bioma_long$nitrogen == nitro )$biomass_per_sp
    
    ### substation between the biomasse of a specie in a triplet anda species in a mono 
    biomass_shift = append(biomass_shift, bioma_long[n,]$biomass_per_sp-
                             biomass_mono)
    
  }
  return(biomass_shift)
}
head(bioma_long)


##### dominante sp ( maybe usefull)

sp_names <- names(cover)[6:17]

cover$dominant_sp = apply(cover[,6:17], 1, get.dominante_sp)[1,]
cover$dominant_strat = apply(cover[,6:17], 1, get.dominante_sp)[2,]


##### biomass shift #####
biomass_per_sp = matrix(ncol = ncol(cover[,6:17]), nrow = nrow(cover) )
for (n in 1:nrow(cover)){
  
  biomass_per_sp[n,] <- biomass.per.sp(bioma[n,], cover[n,])
}

biomass_per_sp = as.data.frame(biomass_per_sp)
names(biomass_per_sp) = names(cover[,6:17])
bioma = cbind(bioma, biomass_per_sp)

bioma_long = bioma %>%
  pivot_longer(cols = Be:To, names_to = "specie_solo", values_to = "biomass_per_sp")
              
str(bioma)
str(interactions)
names(interactions)[1] = "nitrogen"
biomass_shift <- biomass.shift(bioma_long) 
length(biomass_shift)
nrow(bioma_long)
bioma_long = cbind(bioma_long, biomass_shift)

main_data =full_join(bioma_long, interactions, by = c("nitrogen", "species"))
str(main_data)



##### biomass models #####
main_data$feasibility = as.factor(main_data$feasibility)
main_data_control = subset(main_data, main_data$nitrogen == 0 &
                             main_data$type == "triplet")
str(main_data_control)

ggplot(main_data_control)+
  geom_point(aes(Omega, biomass_shift))
### 
str(main_data_control)

mod_biom_shift <- lmer(formula = biomass_shift ~ Omega + feasibility + indirect_effects + log(year_extinction) +(1|specie_solo) , data = main_data_control)
summary(lm(formula = biomass_shift ~ Omega + feasibility + indirect_effects + log(year_extinction) , data = main_data_control))$r.squared

coef_biom_shift = summary(mod_biom_shift)$coefficients
hist(summary(mod_biom_shift)$residuals)

ggplot(main_data_control)+
  geom_point(aes(Omega, biomass_shift,col =plot ))+
  geom_abline(intercept = coef_biom_shift[1,1], slope = coef_biom_shift[2,1],
              linetype =2 , color = "blue")
ggplot(main_data_control)+
  geom_boxplot(aes(feasibility, biomass_shift))+
  geom_abline(intercept = coef_biom_shift[1,1], slope = coef_biom_shift[3,1],
              linetype =2 , color = "blue")

ggplot(main_data_control)+
  geom_point(aes(indirect_effects, biomass_shift))+
  geom_abline(intercept = coef_biom_shift[1,1], slope = coef_biom_shift[4,1],
              linetype =2 , color = "blue")

ggplot(main_data_control)+
  geom_point(aes(year_extinction , biomass_shift))+
  geom_abline(intercept = coef_biom_shift[1,1], slope = coef_biom_shift[5,1],
              linetype =2 , color = "blue")+
  scale_x_log10()
#### the problem is might be that we have 3 biomasse shift value per plot, one for 
# each species. Making a models on those value make no sense because omega, indirect
# effects, ... are values per plot, not per specie. 

b_shift_mean=c()
for(n_plot in 1:48){
  ### get only the bioma_sift of the species of the triplet 
  temp = subset(main_data, main_data$plot == n_plot)
  
  temp1 = subset(temp, temp$specie_solo == temp$species_1[1])
  temp2 = subset(temp, temp$specie_solo == temp$species_2[1])
  temp3 = subset(temp, temp$specie_solo == temp$species_3[1])
  temp = rbind(temp1,rbind(temp2, temp3))
  b_shift_mean = append(b_shift_mean,
                        mean(subset(temp, temp$nitrogen == 0)$biomass_shift,
                             na.rm = TRUE))
  
  b_shift_mean = append(b_shift_mean,
                        mean(subset(temp, temp$nitrogen == 1)$biomass_shift,
                             na.rm = TRUE))
} 

b_shift_mean = data.frame(b_shift_mean, nitrogen = rep(c(0,1), 48*2))
b_shift_mean = cbind(b_shift_mean,
                     species = subset(bioma,
                            as.numeric(as.character(bioma$plot)) <49 )$species)
b_shift_mean$nitrogen =as.factor(b_shift_mean$nitrogen)
str(b_shift_mean)
main_data = full_join(main_data, b_shift_mean , by = c("nitrogen", "species") )

##### same than before but with the mean of biomasse shift
head(main_data_control)
main_data_control = subset(main_data, main_data$nitrogen == 0 &
                             main_data$type == "triplet")
mod_biom_shift2 <- lmer(formula = b_shift_mean ~ Omega + feasibility + indirect_effects + log(year_extinction) +(1|specie_solo), data = main_data_control)
summary(mod_biom_shift2)
coef_biom_shift2 = summary(mod_biom_shift)$coefficients

hist(summary(mod_biom_shift2)$residuals)
hist(main_data_control$b_shift_mean)
ggplot(main_data_control)+
  geom_point(aes(Omega, b_shift_mean))+
  geom_abline(intercept = coef_biom_shift2[1,1], slope = coef_biom_shift2[2,1],
              linetype =2 , color = "blue")
  
ggplot(main_data_control)+
  geom_boxplot(aes(feasibility, b_shift_mean))+
  geom_abline(intercept = coef_biom_shift2[1,1], slope = coef_biom_shift2[3,1],
              linetype =2 , color = "blue")

ggplot(main_data_control)+
  geom_point(aes(indirect_effects, b_shift_mean))+
  geom_abline(intercept = coef_biom_shift2[1,1], slope = coef_biom_shift2[4,1],
              linetype =2 , color = "blue")

ggplot(main_data_control)+
  geom_point(aes(year_extinction , b_shift_mean))+
  geom_abline(intercept = coef_biom_shift2[1,1], slope = coef_biom_shift2[5,1],
              linetype =2 , color = "blue")+
  scale_x_log10()
names_coef = attr(coef_biom_shift2, "dimnames")[[1]][-1]
names_coef[2] = "feasibility"
R2_partitioning = partR2(mod_biom_shift2,names_coef)$R2[1:5, 1:2]

R2_partitioning$sign = ifelse(R2_partitioning$estimate > 0, 
       "positive",
       "negaive")
str(R2_partitioning)
ggplot(R2_partitioning)+
  scale_x_continuous(limits = c(-0.5,0.5))+
  geom_col(aes(estimate, term, fill = sign),width = 0.3)+
  theme_classic()



