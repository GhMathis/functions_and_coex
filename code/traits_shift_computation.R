library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)

library(FactoMineR)
library(vegan)
library(factoextra)
#########################################################
    ########## Shift of functional traits ##########
#########################################################

localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)
getwd()
data_dir = "organized_data"
code_dir ="code"

##### load the data
load(file.path(data_dir,"bioma.RData"))
load(file.path(data_dir,"interactions_triplet.RData"))
load(file.path(data_dir,"cover.RData"))
load(file.path(data_dir,"enzimatic_activity.RData"))
load(file.path(data_dir,"herbivoty_pathogen.RData"))

##### load the functions
source(file.path(code_dir, "traits_shift_functions.R"))

########## shape the data ########## 
cover_long = cover %>%
  pivot_longer(cols = Be:legumes, names_to = "species_cover", values_to = "perc_cover")

str(cover_long)
str(cover)

hist(cover_long$perc_cover)
ggplot(cover_long)+
  geom_point(aes(species_cover,perc_cover))


##### biomass

bioma = bioma[, c("nitrogen", "plot", "biomass", "species_1",
                  "species_2","species_3", "type")]
list_data_bioma = shape.df(bioma,cover,"biomass")
trait_mono_bioma = list_data_bioma[[1]]
trait_triplet_bioma = list_data_bioma[[2]]

##### enzimatic activity

## Beta-gluscosidase
beta_glu_df = subset(enz_act, enz_act$enz =="beta_glu")

nrow(beta_glu_df)
nrow(cover)

list_data_betaG = shape.df(beta_glu_df,cover,"nytro_release")
trait_mono_betaG = list_data_betaG[[1]]
trait_triplet_betaG = list_data_betaG[[2]]
nrow(trait_triplet_betaG)
nrow(trait_mono_betaG)
head(trait_triplet_betaG)
head(trait_mono_betaG)

## Phosphatase
phospha_df = subset(enz_act, enz_act$enz =="phospha")

nrow(phospha_df)
nrow(cover)

list_data_phospha = shape.df(phospha_df,cover,"nytro_release")
trait_mono_phospha = list_data_phospha[[1]]
trait_triplet_phospha = list_data_phospha[[2]]
nrow(trait_triplet_phospha)
nrow(trait_mono_phospha)
head(trait_triplet_phospha)
head(trait_mono_phospha)



########## Biomass shift ########## 


trait_metrics_biomass = commu.traits.metrics(trait_mono = trait_mono_bioma,
                                             trait_triplet = trait_triplet_bioma,
                                             cover = cover,
                                             trait_name_to = "biomass" )

traits_shift_data=c()

str(interactions)
trait_metrics_biomass_control = subset(trait_metrics_biomass,
                                       trait_metrics_biomass$nitrogen == 0)
str(trait_metrics_biomass_control)
traits_shift_data = full_join(trait_metrics_biomass_control,interactions,
                         by = c("species","nitrogen"))
str(traits_shift_data)


########## Enzymatic shift ########## 

## Beta-gluscosidase

trait_metrics_betaG = commu.traits.metrics(trait_mono = trait_mono_betaG,
                                             trait_triplet = trait_triplet_betaG,
                                             cover = cover,
                                             trait_name_to = "beta_glu" )
head(trait_metrics_betaG)

trait_metrics_betaG_control = subset(trait_metrics_betaG,
                                       trait_metrics_betaG$nitrogen == 0)
str(trait_metrics_biomass_control)
traits_shift_betaG = full_join(trait_metrics_betaG_control,interactions,
                              by = c("species","nitrogen"))
traits_shift_data = rbind(traits_shift_data, traits_shift_betaG)
nrow(traits_shift_data)

## Phosphatase

trait_metrics_phospha = commu.traits.metrics(trait_mono = trait_mono_phospha,
                                           trait_triplet = trait_triplet_phospha,
                                           cover = cover,
                                           trait_name_to = "phospha" )
head(trait_metrics_phospha)


trait_metrics_phospha_control = subset(trait_metrics_phospha,
                                     trait_metrics_phospha$nitrogen == 0)
str(trait_metrics_phospha_control)
traits_shift_phospha = full_join(trait_metrics_phospha_control,interactions,
                               by = c("species","nitrogen"))
traits_shift_data = rbind(traits_shift_data, traits_shift_phospha)
nrow(traits_shift_data)
head(traits_shift_data)



save(traits_shift_data,file = "organized_data/traits_shift_data.RData")

########## Multivariate visualization ##########

biomass_shift_pca<-PCA(subset(traits_shift_data[,c(3:6, 11:18)],
                            traits_shift_data$trait == "biomass" ),quanti.sup = 5:12)
get_eigenvalue(biomass_shift_pca)
fviz_eig(biomass_shift_pca, adlabels=TRUE)
fviz_pca_var(biomass_shift_pca, col.var = "black")
fviz_pca_ind(biomass_shift_pca, col.var = "black")
fviz_pca_biplot(biomass_shift_pca, col.var = "black")
str(biomass_shift_pca)
plot.PCAbiomass_shift_pca