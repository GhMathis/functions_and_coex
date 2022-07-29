library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
###################################################
    ########## ENZYMATIQUE ANNALYSYS ##########
###################################################

localDir="~/Fac/Cesure2/Plant species coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "organized_data"

load(file.path(data_dir,"enzimatic_activity.RData"))
load(file.path(data_dir,"interactions_triplet.RData"))
load(file.path(data_dir,"coex.RData"))
str(enz_act)
str(interactions)


##### some visuals verification of the wet and frost soils #####
### cali_group = soil samples incubate during the same day
enz_act2 <- read.delim2("oster_data_2022/enz_act.txt", header = TRUE)
str(enz_act2)
temp = subset(enz_act2, sample == "a" & plot != 0)
str(temp)
enz_act$cali_group = as.factor(temp$cali_group)
enz_act$comments
ggplot(subset(enz_act, comments != ""),
       aes(plot, nytro_release, col = enz, label = comments))+
  geom_text()
mean(subset(enz_act, comments != "" & type == "triplet" & enz == "beta_glu")$nytro_release)
mean(subset(enz_act, comments == "" & type == "triplet" & enz == "beta_glu")$nytro_release)
mean(subset(enz_act, comments != "" & type == "triplet" & enz == "phospha")$nytro_release)
mean(subset(enz_act, comments == "" & type == "triplet" & enz == "phospha")$nytro_release)

##### linear model between control an nitorgen treatment for both enzyme #####
## beta
beta_data <- subset(enz_act, enz == "beta_glu")
model_beta <- lm(beta_data$nytro_release~beta_data$nitrogen)
summary(model_beta)
anova(model_beta)
hist(subset(enz_act$nytro_release, enz_act$enz == "beta_glu"))

par(mfrow = c(2,2))
plot(model_beta, las = 1) ### normal distib OK
par(mfrow = c(1,1))

ggplot(beta_data,aes(nitrogen,nytro_release))+
  geom_boxplot()+
  geom_abline(intercept= 0.70700, slope =-0.01201,
              col = "blue", linetype =2, cex =1.5)
## phospha
phospha_data <- subset(enz_act, enz == "phospha")
model_phospha <- lm(phospha_data$nytro_release~phospha_data$nitrogen)
str(summary(model_phospha))
anova(model_phospha)

par(mfrow = c(2,2))
plot(model_phospha, las = 1) ### normal distib OK
par(mfrow = c(1,1))

ggplot(phospha_data,aes(nitrogen,nytro_release))+
  geom_boxplot()+
  geom_abline(intercept= 1.96564, slope =-0.02553,
              col = "blue", linetype =2, cex =1.5)
### no diff in enzymatique activity between N and C 


##### linear model between triplet an monoculutre for both enzyme #####
## beta
ggplot(enz_act, aes(type,nytro_release,col=nitrogen))+
  geom_boxplot()+
  facet_wrap(~enz)

model_beta_cult <- lm(nytro_release~type,data =beta_data)
summary(model_beta_cult)
anova(model_beta_cult)

par(mfrow = c(2,2))
plot(model_beta_cult, las = 1) ### normal distib OK
par(mfrow = c(1,1))

coef = model_beta_cult$coefficients
str(coef)

## phospha
phospha_data <- subset(enz_act, enz == "phospha")
model_phospha_cult <- lm(nytro_release~type,data=phospha_data)
summary(model_phospha)
anova(model_phospha)

par(mfrow = c(2,2))
plot(model_phospha, las = 1) ### normal distib OK
par(mfrow = c(1,1))

##### test the impact of incubation date on enz. act.#####

## beta
ggplot(enz_act)+
  geom_boxplot(aes(cali_group,nytro_release, col = enz))+
  facet_grid(~nitrogen)
str(beta_data)
model_beta2 <- lm(nytro_release~cali_group+0, data = beta_data)
summary(model_beta2)
anova(model_beta2)

par(mfrow = c(2,2))
plot(model_beta2, las = 1) ### normal distib OK
par(mfrow = c(1,1))

## phosph
model_phospha2 <- lm(nytro_release~cali_group+0, data = phospha_data)
summary(model_phospha2)
anova(model_phospha2)

#### De grosses diferences entre les groupes de calibration.
##
par(mfrow = c(2,2))
plot(model_phospha2, las = 1) ### normal distib OK
par(mfrow = c(1,1))


######## Enzyme and triplet coexistence ########
tail(interactions)
tail(enz_act)
temp = c()
names(interactions)[1]="nitrogen"
interactions$nitrogen = as.factor(rep(0,nrow(interactions)))
coex$species = apply( coex[,3:5], 1, 
        function(x) paste(x[!is.na(x)], collapse = "_"))


enz_act_wider = enz_act[,-c(5,6,8)]
enz_act_wider = enz_act_wider %>%
  pivot_wider(names_from = enz, values_from = nytro_release)
enz_act_wider$species = coex$species
enz_triplet = subset(enz_act_wider, type == "triplet")

str(enz_act_wider)

main_data = full_join(enz_triplet,interactions, by = c("species"))

main_data = main_data %>%
  pivot_longer(c(beta_glu,phospha),
               names_to = "enz", values_to = "nytro_release" )
tail(main_data)
head(main_data)
str(main_data)
nrow(main_data)

ggplot(main_data)+
  facet_wrap(~enz)+
  geom_point(aes(nytro_release, Omega, col= species))

##### 

shape_df <- shape.df(enz_act, cover, "nytro_release")
trait_mono_beta = subset(shape_df[[1]], shape_df[[1]]$enz == "beta_glu")
trait_triplet_beta =subset(shape_df[[2]], shape_df[[2]]$enz == "beta_glu")

commu.traits.metrics(trait_mono = trait_mono_beta, trait_triplet = trait_triplet_beta,
                     cover = cover, trait_name = "nytro_release" )