library(ggplot2)
library(readxl)
library(writexl)

#### load and shape data ####
setwd("~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/data")

enz_act <- read.delim2("enzimatic_activity.txt", header= TRUE, dec = ".",sep=" ")


str(enz_act)


#calculate the dry weight ratio(dry mass/fresh mass), fresh mass was weighted
#precisely around 5g +-0.005g
enz_act$dry_weight_ratio = enz_act$dry_weight/rep(5,nrow(enz_act))

#### find the calibrations regressions ####

### only calibration data
lm_cali = data.frame()
data_cali = subset(enz_act, sample == "cali")
str(data_cali)
for (n_cali in 1:6){
  temp = subset(data_cali, cali_group == n_cali)
  coefs = as.data.frame(lm(Abs~p_nytro_C, data = temp)$coefficients)
  lm_cali = rbind(lm_cali,data.frame(intercept= coefs[1,], slope = coefs[2,],
                                     cali_group = n_cali))
} 
str(lm_cali)

### calculation of the concentrations given linear models of calibrations
data_enz = subset(enz_act, sample != "cali")

for (n_abs in 1:nrow(data_enz)){
  n_cali = enz_act$cali_group[n_abs]
  data_enz$p_nytro_C[n_abs] = ((data_enz$Abs[n_abs] - lm_cali[n_cali,1]) /
                                 (lm_cali[n_cali,2] * data_enz$dry_weight_ratio[n_abs]))
  
  
}
tail(data_enz)
ggplot() +
  geom_point(data = data_enz,
             aes(p_nytro_C, Abs,col=enz,shape=sample),
             cex = 1.5,color = "black")+
  geom_point(data = data_enz,aes(p_nytro_C, Abs,col=enz,shape=sample),cex = 1)+
  facet_wrap(~cali_group)+
  geom_smooth(data =subset(enz_act, sample == "cali"),
              aes(p_nytro_C,Abs ),method = lm,se = F,linetype=2, size=1)+
  scale_y_continuous(limits=c(-0.1,1.1))+
  scale_x_continuous(limits=c(-0.1,350))+
  xlab("enz (in mmol.g-1.h-1)")+
  theme(axis.text.x = element_text( size=25),
        axis.text.y = element_text( size=25),
        axis.title=element_text(size=30))
ggplot()  +
  geom_boxplot(data = data_enz,aes(sample, p_nytro_C))+
  facet_wrap(~enz)
ggplot()+
  facet_wrap(~enz)+
  geom_boxplot(data = data_enz,aes(date_time, p_nytro_C))
head(data_enz)

### sample A - sample B 
nytro_release = (subset(data_enz,sample == "a")$p_nytro_C -
                subset(data_enz,sample == "b")$p_nytro_C) * 1/139
data_enz$nytro_release = nytro_release
enz_act = cbind(nytro_release, subset(data_enz,sample == "b",
                                select = c(enz, nitrogen, plot,comments,date_time)))

str(enz_act)
str(data_enz)
enz_act$nytro_release = NA
data_enz = rbind(data_enz, subset(enz_act, enz_act$sample == "cali"))

write.table(data_enz, "enzimatic_activity.txt", dec = ".", col.names = TRUE)
