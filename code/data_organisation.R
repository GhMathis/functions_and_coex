install.packages("dplyr")
install.packages("tidyr")
for (i in 1){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(writexl)
}

#### load and shape data ####
setwd("~/Fac/Cesure2/Plant species coexistence/annalisys/oster_data_2022")
enz_act <- read.delim2("enz_act.txt", header = TRUE)
soil_weight <- read.delim2("soil_weight.txt", header = TRUE)
soil_timing <- read.delim2("soil_timing.txt", header = TRUE)
par_height <- read.delim2("PAR_and_height6.txt",header = TRUE)
herbi <- read.delim2("herbi_patho_ostermundigen_2022.txt", header = TRUE)
coex <-read.delim2("coexistence_2022.txt", header =TRUE)
###### Enzymatic data #######
head(soil_weight)
str(soil_weight)
head(enz_act)

### format, shape and merge the data properly
## enz_act
enz_act = enz_act[,1:7]
enz_act$sample[enz_act$sample=="na"]="cali"
enz_act[enz_act=="na"]=NA
enz_act=enz_act[!is.na(enz_act$plot),]
enz_act$Abs = as.numeric(enz_act$Abs)
enz_act$nitrogen = as.factor(enz_act$nitrogen)
enz_act$plot = as.factor(enz_act$plot)
enz_act$cali_group = as.factor(enz_act$cali_group)
names(enz_act)[6]="p_nytro_C"
enz_act$p_nytro_C = enz_act$p_nytro_C

##soil_weight
names(soil_weight)[1:2]=c("plot","nitrogen")
soil_weight=soil_weight[!is.na(soil_weight$plot),]
soil_weight$plot = as.factor(soil_weight$plot)
soil_weight$nitrogen = as.factor(soil_weight$nitrogen)
soil_weight$dry_weight = as.numeric(soil_weight$dry_weight)

## soil_timing
head(soil_timing)
str(soil_timing)
names(soil_timing)[1:3] = c("plot","nitrogen","date_time")
soil_timing$enz = rep("phospha",120)
soil_timing$plot = as.factor(soil_timing$plot)
soil_timing$nitrogen[soil_timing$nitrogen == "C"] = 0
soil_timing$nitrogen[soil_timing$nitrogen == "N"] = 1
soil_timing$nitrogen = as.factor(soil_timing$nitrogen)

#calculate the dry weight ratio(dry mass/fresh mass), fresh mass was weighted
#precisely around 5g +-0.005g
soil_weight$dry_weight_ratio = soil_weight$dry_weight/rep(5,nrow(soil_weight))

# merge data
enz_act = full_join(enz_act,soil_weight, by = c("plot", "nitrogen"))
enz_act = full_join(enz_act,soil_timing, by = c("plot", "nitrogen","enz"))

#remove useless col
enz_act = enz_act[, c(1:7, 10:13)]

head(enz_act)
str(enz_act)
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
data = subset(enz_act, sample != "cali")
for (n_abs in 1:nrow(data)){
    n_cali = enz_act$cali_group[n_abs]
    data$p_nytro_C[n_abs] = ((data$Abs[n_abs] - lm_cali[n_cali,1]) /
      (lm_cali[n_cali,2] * data$dry_weight_ratio[n_abs]))
      
  
}
tail(data[480,])
ggplot() +
  geom_point(data = data,
             aes(p_nytro_C, Abs,col=enz,shape=sample),
             cex = 1.5,color = "black")+
  geom_point(data = data,aes(p_nytro_C, Abs,col=enz,shape=sample),cex = 1)+
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
  geom_boxplot(data = data,aes(sample, p_nytro_C))+
  facet_wrap(~enz)
ggplot()+
  facet_wrap(~enz)+
  geom_boxplot(data = data,aes(date_time, p_nytro_C))
head(data)

main_data =data.frame(nytro_release = (subset(data,sample == "a")$p_nytro_C -
               subset(data,sample == "b")$p_nytro_C) * 1/139)

main_data = cbind(main_data, subset(data,sample == "b",
                                    select = c(enz, nitrogen, plot,comments)))
subset(main_data, nitrogen == 1 & plot == 22)$comment = c("wet","wet")
subset(main_data, nitrogen == 1 & plot == 23)$comment =  c("wet","wet")
subset(main_data, nitrogen == 1 & plot == 27)$comment =  c("wet","wet")
subset(main_data, nitrogen == 1 & plot == 34)$comment =  c("wet","wet") 
subset(main_data, nitrogen == 1 & plot == 31)$comment =  c("wet","wet") 
subset(main_data, nitrogen == 1 & plot == 11)$comment =  c("wet","wet") 
subset(main_data, nitrogen == 0 & plot == 6)$comment  =  c("wet","wet") 
subset(main_data, nitrogen == 0 & plot == 58)$comment =  c("wet","wet")
head(main_data)
str(main_data)

ggplot(main_data)+
  geom_boxplot(aes(nitrogen,nytro_release))+
  facet_wrap(~enz)
##### PAR and height data #####

str(par_height)
nrow(par_height)
names(par_height)[c(1,11)] = c("plot", "nitrogen")
par_height$date_char = par_height$date 
par_height$date = as.Date(par_height$date, format="%m.%d.%y")

for (i in 5:10){
  par_height[,i] = as.numeric(par_height[,i])
}
for (i in c(1,11,15) ){
  par_height[,i] = as.factor(par_height[,i])
}


par_height$PAR_inf[par_height$PAR_inf>21] = round(par_height$PAR_inf[par_height$PAR_inf>21])
par_height$mean_height = rowMeans(par_height[,5:8])
par_height$PAR_diff = par_height$PAR_sup- par_height$PAR_inf

str(par_height)

par_height=par_height%>%
  pivot_longer(c(species_1, species_2, species_3), names_to = "attribute",
               values_to = "species")
par_height=par_height[c(par_height$species != "-"),]
str(par_height)

ggplot(par_height)+
  facet_wrap(~nitrogen)+
  geom_boxplot(aes(date_char, mean_height))
  geom_smooth(aes(date, mean_height), se=FALSE, span = 0.9)
  
##### herbi data#####

head(herbi)

str(herbi)
nrow(herbi)
for(i in c(2,3)){
  herbi[,i] = as.factor(herbi[,i])
}
for(i in 5:20){
  herbi[,i] = as.numeric(herbi[,i])
}

ggplot(herbi)+
  geom_boxplot(aes(species,perc_dmg))

##### recap data #####
str(main_data)
main_data = main_data%>%
  pivot_wider(names_from = enz, values_from = nytro_release)

nrow(main_data)
ggplot(main_data)+
  geom_point(aes(plot, nytro_release, col = enz))
head(par_height)
nrow(par_height)
head(herbi)
nrow(herbi)
head(coex)
nrow(coex)
