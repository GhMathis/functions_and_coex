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

localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "organized_data"

load(file.path(data_dir,"bioma.RData"))
load(file.path(data_dir,"interactions_triplet.RData"))
load(file.path(data_dir,"cover.RData"))
head(bioma)

##### shape data #####

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
head(bioma)

##### functions #####

trait.per.sp <- function(trait_df, cover, trait_name){
  #### function to find the biomass per species with the biomass data
  # and the cover data.
  names(trait_df)[names(trait_df) == trait_name] = "trait"

  # cover weighted biomass calcul 
  trait_weighted = (cover[,c("Be", "Cb", "Dc", "Fr", "Hl","Lp", "Pg",
                             "Pm", "Pt", "Ra", "Sp", "To")]/cover$total.cover) * trait_df$trait
  return(as.matrix(trait_weighted))
}


commu.traits.metrics <- function (trait_mono, trait_triplet, cover, trait_name,
                                  n_sp = 3, n_plot = 96 ){
  #### Function to get 4 trait metrics per plot : The sown trait value, delta
  # abundance shift, delta intraspecific shift and CWM

  ### trait_mono : wide df of the trait values for species in mono culture
  # it need to have cols of : nitrogen, species_1 and the trait values of the specie.
  
  ### trait_triplet : wide df of the trait values for species in triplet
  # it need to have cols of : nitrogen, species_1, species_2, species_3
  # plot and the 12 col of each species trait values ("Be", "Cb", "Dc", "Fr", "Hl",
  # "Lp", "Pg", "Pm", "Pt", "Ra", "Sp", "To")
  
  ### cover : wide df  with cover of each sp and the total cover. Cols :
  # sp("Be", "Cb", "Dc", "Fr", "Hl","Lp", "Pg", "Pm", "Pt", "Ra", "Sp", "To") 
  # and total.cover
  
  ### n : number of species in the plot (always 3 in our case)
  ### n_plot : number of plot (always 96 in our case)
  ### trait_name : the name of the trait col in the df "trait_mono".
  
  trait_metrics = data.frame(nitrogen = trait_triplet$nitrogen, plot = trait_triplet$plot)
  
  ##### relative abundance
  relative_ab =c()
  ### cover_sp/cover_tot
  for (n in 1:nrow(cover)){
    relative_ab = rbind(relative_ab, cover[n,c("Be", "Cb", "Dc", "Fr", "Hl",
                                               "Lp", "Pg", "Pm", "Pt", "Ra", "Sp", "To")]/cover$total.cover[n])
  }
  relative_ab_triplet = as.matrix(relative_ab[1:96,])

  control_trait_mat = matrix(data= NA, nrow = n_plot, ncol = 12)
  colnames(control_trait_mat) = colnames(relative_ab_triplet)
  
  #### sown trait per plot
  sown_trait = c()
  
  ### sum(sp_trait_mono)/n_sp_plot
  for (n in 1:nrow(trait_triplet)){
    sp_trait = c()
    
    sp_trait = append(sp_trait, subset(trait_mono,
                                       trait_mono$species_1 == trait_triplet$species_1[n] &
                                         trait_mono$nitrogen == trait_triplet$nitrogen[n])$trait)
    
    sp_trait = append(sp_trait, subset(trait_mono,
                                        trait_mono$species_1 == trait_triplet$species_2[n] &
                                          trait_mono$nitrogen == trait_triplet$nitrogen[n])$trait)
    
    sp_trait = append(sp_trait, subset(trait_mono,
                                        trait_mono$species_1 == trait_triplet$species_3[n] &
                                          trait_mono$nitrogen == trait_triplet$nitrogen[n])$trait)
    
    sp_names = trait_triplet[n,"species_1"]
    sp_names = append(sp_names, trait_triplet[n,"species_2"])
    sp_names = append(sp_names, trait_triplet[n,"species_3"])
    
    control_trait_mat[n,sp_names] = sp_trait
    
    sown_trait = append(sown_trait, sum(sp_trait, na.rm = T))/n_sp
  }
  
  trait_metrics$sown_trait = sown_trait
 
  
  ##### delta abondance shift per plot 
  ### sum(Ab_rela*sp_trait_mono) - sown_trait
  
  mono_weight_mean= relative_ab_triplet * control_trait_mat 
  delta_ab_shift = rowSums(mono_weight_mean, na.rm = T) - sown_trait
  trait_metrics$delta_ab_shift = delta_ab_shift
  
  
  ##### delta intraspecific shift per plot
  ### sum(sp_trait_triplet)/n_sp_plot - sown_trait
  
  trait_values = as.matrix(trait_triplet[, c("Be", "Cb", "Dc", "Fr", "Hl",
                                  "Lp", "Pg", "Pm", "Pt", "Ra", "Sp", "To")])
  delta_intra_shift = rowSums(trait_values, na.rm = T) / n_sp - sown_trait
  trait_metrics$delta_intra_shift = delta_intra_shift
  
  
  ##### Comunity weighted mean
  ###sum(Ab_rela * sp_trait_triplet)
  CWM = rowSums(relative_ab_triplet * trait_values, na.rm = T)
  trait_metrics$CWM = CWM
  trait_metrics$species = trait_triplet$species
  
  trait_metrics$trait =  rep(trait_name, nrow(trait_metrics))
  return(trait_metrics)
}


shape.df <- function(df,cover,trait_name){
  
  
  if (!any(names(df) == "species_1" )){
    df = cbind(df, cover[, c("species_1", "species_2", "species_3")])
  }
  
  df$species = apply( df[,c("species_1", "species_2", "species_3")], 1, 
                      function(x) paste(x[!is.na(x)], collapse = "_"))
  trait_per_sp = trait.per.sp(df, cover, trait_name )
  df = cbind(df, trait_per_sp)
  
  trait_triplet = subset(df, df$type == "triplet")
  trait_mono = subset(df, df$type == "mono")
  
  ### get the trait value of the species in solo (without the alien species) 
  n=1
  trait = c()
  for (n in 1:nrow(trait_mono)){
    trait = append(trait, trait_mono[n, c(trait_mono$species_1[n]) ])
  }
  trait_mono$trait = trait


  trait_mono[,c("Be", "Cb", "Dc", "Fr", "Hl", "Lp", "Pg",
                "Pm", "Pt", "Ra", "Sp", "To", "species_2", "species_3")] <- NULL  
  return(list(trait_mono = trait_mono, trait_triplet = trait_triplet))
}
##### shape data and find the trait metrics of each plot for the biomass #####

nrow(bioma)
bioma = bioma[, 1:7]
list_data_bioma = shape.df(bioma,cover,"biomass")
trait_mono_bioma = list_data_bioma[[1]]
trait_triplet_bioma = list_data_bioma[[2]]

trait_metrics_biomass = commu.traits.metrics(trait_mono = trait_mono_bioma,
                                             trait_triplet = trait_triplet_bioma,
                                             cover = cover,
                                             trait_name = "biomass" )

biomass_data=c()

str(interactions)
trait_metrics_biomass_control = subset(trait_metrics_biomass,
                                       trait_metrics_biomass$nitrogen == 0)
str(trait_metrics_biomass_control)
biomass_data = full_join(trait_metrics_biomass_control,interactions,
                         by = c("species","nitrogen"))
str(biomass_data)
biomass_data$feasibility = as.factor(biomass_data$feasibility)
##### models CWM #####
hist(biomass_data$CWM,breaks =48 )
hist(biomass_data$delta_ab_shift,breaks =48 )
hist(biomass_data$delta_intra_shift,breaks =48 )
mod_CWM <- lm(formula = CWM ~ Omega + feasibility + indirect_effects + log(year_extinction), data = biomass_data)
summary(mod_CWM)
coef_CWM = summary(mod_CWM)$coefficients

hist(summary(mod_CWM)$residuals, breaks = 48)
hist(biomass_data$CWM)
ggplot(biomass_data)+
  geom_point(aes(Omega, CWM))+
  geom_abline(intercept = coef_CWM[1,1], slope = coef_CWM[2,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_boxplot(aes(feasibility, CWM))+
  geom_abline(intercept = coef_CWM[1,1], slope = coef_CWM[3,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_point(aes(indirect_effects, CWM))+
  geom_abline(intercept = coef_CWM[1,1], slope = coef_CWM[4,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_point(aes(year_extinction , CWM))+
  geom_abline(intercept = coef_CWM[1,1], slope = coef_CWM[5,1],
              linetype =2 , color = "blue")+
  scale_x_log10()

################### for later ############################
names_coef = attr(coef_CWM, "dimnames")[[1]][-1]
names_coef[2] = "feasibility"
R2_partitioning = partR2(mod_CWM,names_coef)$R2[1:5, 1:2]

R2_partitioning$sign = ifelse(R2_partitioning$estimate > 0, 
                              "positive",
                              "negaive")
str(mod_CWM)
str(R2_partitioning)
ggplot(R2_partitioning)+
  scale_x_continuous(limits = c(-0.5,0.5))+
  geom_col(aes(estimate, term, fill = sign),width = 0.3)+
  theme_classic()
##########################################################


##### delta ab shift #####

mod_delta_ab_shift <- lm(formula = delta_ab_shift ~ Omega + feasibility + indirect_effects + log(year_extinction), data = biomass_data)
summary(mod_delta_ab_shift)
coef_delta_ab_shift = summary(mod_delta_ab_shift)$coefficients

hist(summary(mod_delta_ab_shift)$residuals, breaks = 48)
hist(biomass_data$delta_ab_shift)
ggplot(biomass_data)+
  geom_point(aes(Omega, delta_ab_shift))+
  geom_abline(intercept = coef_delta_ab_shift[1,1], slope = coef_delta_ab_shift[2,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_boxplot(aes(feasibility, delta_ab_shift))+
  geom_point(aes(feasibility, delta_ab_shift))+
  geom_abline(intercept = coef_delta_ab_shift[1,1], slope = coef_delta_ab_shift[3,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_point(aes(indirect_effects, delta_ab_shift))+
  geom_abline(intercept = coef_delta_ab_shift[1,1], slope = coef_delta_ab_shift[4,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_point(aes(year_extinction , delta_ab_shift))+
  geom_abline(intercept = coef_delta_ab_shift[1,1], slope = coef_delta_ab_shift[5,1],
              linetype =2 , color = "blue")+
  scale_x_log10()

##### delta intra sp shift #####

mod_delta_intra_shift <- lm(formula = delta_intra_shift ~ Omega + feasibility + indirect_effects + log(year_extinction), data = biomass_data)
summary(mod_delta_intra_shift)
coef_delta_intra_shift = summary(mod_delta_intra_shift)$coefficients

hist(summary(mod_delta_intra_shift)$residuals, breaks = 48)
hist(biomass_data$delta_intra_shift)
ggplot(biomass_data)+
  geom_point(aes(Omega, delta_intra_shift))+
  geom_abline(intercept = coef_delta_intra_shift[1,1], slope = coef_delta_intra_shift[2,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_boxplot(aes(feasibility, delta_intra_shift))+
  geom_point(aes(feasibility, delta_intra_shift))+
  geom_abline(intercept = coef_delta_intra_shift[1,1], slope = coef_delta_intra_shift[3,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_point(aes(indirect_effects, delta_intra_shift))+
  geom_abline(intercept = coef_delta_intra_shift[1,1], slope = coef_delta_intra_shift[4,1],
              linetype =2 , color = "blue")

ggplot(biomass_data)+
  geom_point(aes(year_extinction , delta_intra_shift))+
  geom_abline(intercept = coef_delta_intra_shift[1,1], slope = coef_delta_intra_shift[5,1],
              linetype =2 , color = "blue")+
  scale_x_log10()
