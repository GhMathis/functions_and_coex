library(ggplot2)
library(readxl)
library(writexl)

#### load and shape data ####
setwd("~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/oster_data_2022")
enz_act <- read.delim2("enz_act.txt",header= TRUE)
enz_act = enz_act[,1:7]
tail(enz_act)

str(enz_act)
enz_act$Abs = as.numeric(enz_act$Abs)
enz_act$nitrogen = as.integer(enz_act$nitrogen)
enz_act$plot = as.factor(enz_act$plot)
enz_act$cali_group = as.factor(enz_act$cali_group)

#### find the calibrations regressions ####

lm_cali = data.frame()
data_cali = subset(enz_act, sample == "cali")
str(data_cali)
for (n_cali in 1:6){
  print(n_cali)
  temp = subset(data_cali, cali_group == n_cali)
  coefs = as.data.frame(lm(cali_cons~Abs, data = temp)$coefficients)
  lm_cali = rbind(lm_cali,data.frame(intercept= coefs[1,], slope = coefs[2,],
                              cali_group = n_cali))
} 

data = subset(enz_act, sample != "cali")

for (n_abs in 1:nrow(data)){
    n_cali=enz_act$cali_group[n_abs]
    data$cali_cons[n_abs]=lm_cali[n_cali,2]*enz_act$Abs[n_abs]+lm_cali[n_cali,1]
  
}

str(lm_cali)
plot(lm_cali[[1]])
str(data)
ggplot() +
  geom_point(data = data,
             aes(Abs, cali_cons,col=enz,shape=sample))+
  facet_wrap(~cali_group)+
  geom_point(data =subset(enz_act, sample == "cali"),
              aes(Abs, cali_cons),color ="blue")+
  geom_smooth(data =subset(enz_act, sample == "cali"),
              aes(Abs, cali_cons),method = lm,se = F)
