for (i in 1){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(writexl)
}

#### load and shape data ####
setwd("~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/oster_data_2022")
enz_act <- read.delim2("enz_act.txt", header = TRUE)
soil_weight <- read.delim2("soil_weight.txt", header = TRUE)
soil_timing <- read.delim2("soil_timing.txt", header = TRUE)
par_height <- read.delim2("PAR_and_height6.txt",header = TRUE)
herbi <- read.delim2("herbi_patho_ostermundigen_2022.txt", header = TRUE)
cover <-read.delim2("cover_2022.txt", header =TRUE)
bioma <- read.delim2("biomasse_OM_2022.txt", header =TRUE)
###### Enzymatic data #######
head(soil_weight)
str(soil_weight)
head(enz_act)
test = 
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
enz_act = enz_act[order(enz_act$enz, enz_act$plot),]

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
enz_act_brute_data = enz_act[, c(1:7, 10:13)]
enz_act_brute_data = enz_act

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

enz_act =data.frame(nytro_release = (subset(data_enz,sample == "a")$p_nytro_C -
                                       subset(data_enz,sample == "b")$p_nytro_C) * 1/139)

enz_act = cbind(enz_act, subset(data_enz,sample == "b",
                                select = c(enz, nitrogen, plot,comments,date_time)))
enz_act

enz_act[c(enz_act$nitrogen == 1 & enz_act$plot == 22),]$comments = c("wet","wet")
enz_act[c(enz_act$nitrogen == 1 & enz_act$plot == 23),]$comments = c("wet","wet")
enz_act[c(enz_act$nitrogen == 1 & enz_act$plot == 27),]$comments = c("wet","wet")
enz_act[c(enz_act$nitrogen == 1 & enz_act$plot == 34),]$comments = c("wet","wet")
enz_act[c(enz_act$nitrogen == 1 & enz_act$plot == 31),]$comments = c("wet","wet")
enz_act[c(enz_act$nitrogen == 1 & enz_act$plot == 11),]$comments = c("wet","wet")
enz_act[c(enz_act$nitrogen == 0 & enz_act$plot == 6),]$comments = c("wet","wet")
enz_act[c(enz_act$nitrogen == 0 & enz_act$plot == 58),]$comments = c("wet","wet")

enz_act$type = ifelse(as.numeric(as.character(enz_act$plot)) <= 48 ,
                      "triplet", "mono")

head(enz_act)
str(enz_act)

ggplot(enz_act)+
  geom_boxplot(aes(nitrogen,nytro_release))+
  facet_wrap(~enz)
##### PAR and height data #####
par_height= par_height[,-c(16,17)]
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

par_height$species_2[c(par_height$species_2 == "-")] = NA
par_height$species_3[c(par_height$species_3 == "-")] = NA
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
##### cover data #####

cover$species.2[cover$species.2 == "-"] = NA
cover$species.3[cover$species.3 == "-"] = NA
names(cover)[1 : 5]= c("plot", "nitrogen", "species_1", "species_2", "species_3")

##### biomasse data #####
str(bioma)
bioma$biomasse = as.numeric(bioma$biomasse)
bioma$plot = as.factor(bioma$plot)
bioma$nitrogen = as.factor(bioma$nitrogen)
bioma$type = ifelse(as.numeric(as.character(bioma$plot)) <= 48 ,
                      "triplet", "mono")
bioma <- bioma  %>% mutate_all(na_if,"")
##### recap data #####
str(enz_act)
nrow(enz_act)
ggplot(enz_act)+
  geom_point(aes(plot, nytro_release, col = enz))
head(par_height)
nrow(par_height)
head(herbi)
nrow(herbi)
head(cover)
nrow(cover)

tail(bioma)
nrow(bioma)

##### Rdata #####

setwd("~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex")

save(enz_act,file = "organized_data/enzimatic_activity.RData")
save(par_height,file = "organized_data/par_and_height.RData")
save(herbi,file = "organized_data/herbivoty_pathogen.RData")
save(cover,file = "organized_data/cover.RData")
save(bioma,file = "organized_data/bioma.RData")
save(interactions, file = "organized_data/interactions_triplet.RData")
nrow(enz_act)
