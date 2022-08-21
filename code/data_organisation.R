for (i in 1){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(writexl)
}

#### load and shape data ####
localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "oster_data_2022"


enz_act <- read.delim2(file.path(data_dir,"enz_act.txt"), header = TRUE, dec = ".")
soil_weight <- read.delim2(file.path(data_dir,"soil_weight.txt"), header = TRUE, dec = ".")
soil_timing <- read.delim2(file.path(data_dir,"soil_timing.txt"), header = TRUE, dec = ".")
par_height <- read.delim2(file.path(data_dir,"PAR_and_height6.txt"),header = TRUE, dec = ".")
herbi <- read.delim2(file.path(data_dir,"herbi_patho_ostermundigen_2022.txt"), header = TRUE, dec = ".")
cover <-read.delim2(file.path(data_dir,"cover_2022.txt"), header =TRUE, dec = ".")
bioma <- read.delim2(file.path(data_dir,"biomasse_OM_2022.txt"), header =TRUE, dec = ".")
###### Enzymatic data #######
head(soil_weight)
str(soil_weight)
head(enz_act)
tail(enz_act)
str(enz_act)
### format, shape and merge the data properly
## enz_act

enz_act[enz_act=="na"]=NA
enz_act=enz_act[!is.na(enz_act$plot),]

enz_act$Abs = as.numeric(enz_act$Abs)

names(enz_act)[6]="p_nytro_C"

enz_act$nitrogen = as.character(enz_act$nitrogen)
enz_act$nitrogen[enz_act$nitrogen == "1"] = "nitrogen"
enz_act$nitrogen[enz_act$nitrogen == "0"] = "control"
names(enz_act)[names(enz_act) == "plot" | names(enz_act) == "nitrogen"] = c("treatment", "plot_ID")

enz_act = enz_act[order(enz_act$enz, enz_act$plot),]

##soil_weight
soil_weight$treatement = as.character(soil_weight$treatement)
soil_weight$treatement[soil_weight$treatement == "1"] = "nitrogen"
soil_weight$treatement[soil_weight$treatement == "0"] = "control"


str(soil_weight)
soil_weight=soil_weight[!is.na(soil_weight$plot),]
soil_weight$plot = as.factor(soil_weight$plot)
soil_weight$nitrogen = as.factor(soil_weight$nitrogen)
soil_weight$dry_weight = as.numeric(soil_weight$dry_weight)

## soil_timing
head(soil_timing)
str(soil_timing)
names(soil_timing)[2]= "treatment"
soil_timing$treatment[soil_timing$treatment == "N"] = "nitrogen"
soil_timing$treatment[soil_timing$treatment == "C"] = "control"
soil_timing$enz = rep("phospha",120)

names(soil_weight)[2]="treatment"


# merge data
enz_act = full_join(enz_act,soil_weight, by = c("plot_ID", "treatment"))
enz_act = full_join(enz_act,soil_timing, by = c("plot_ID", "treatment","enz"))
head(enz_act)
str(enz_act)
enz_act$comments
enz_act[c(enz_act$treatment == "nitrogen" & enz_act$plot_ID == 22),]$comments = c("wet","wet")
enz_act[c(enz_act$treatment == "nitrogen" & enz_act$plot_ID  == 23),]$comments = c("wet","wet")
enz_act[c(enz_act$treatment == "nitrogen" & enz_act$plot_ID  == 27),]$comments = c("wet","wet")
enz_act[c(enz_act$treatment == "nitrogen" & enz_act$plot_ID  == 34),]$comments = c("wet","wet")
enz_act[c(enz_act$treatment == "nitrogen" & enz_act$plot_ID  == 31),]$comments = c("wet","wet")
enz_act[c(enz_act$treatment == "nitrogen" & enz_act$plot_ID  == 11),]$comments = c("wet","wet")
enz_act[c(enz_act$treatment == "control" & enz_act$plot_ID  == 6),]$comments = c("wet","wet")
enz_act[c(enz_act$treatment == "control" & enz_act$plot_ID == 58),]$comments = c("wet","wet")

enz_act$type = ifelse(as.numeric(as.character(enz_act$plot_ID)) <= 48 ,
                      "triplet", "mono")

head(enz_act)
str(enz_act)

ggplot(enz_act[enz_act$sample!="cali",])+
  geom_boxplot(aes(treatment,Abs))+
  geom_point(aes(treatment,Abs))+
  facet_wrap(~enz)

##### PAR and height data #####
par_height= par_height[,-c(16,17)]
str(par_height)
tail(par_height)
nrow(par_height)
names(par_height)[names(par_height) == "treatement"] =  "treatment"
par_height$treatment = as.character(par_height$treatment)
par_height$treatment[par_height$treatment == "1"] = "nitrogen"
par_height$treatment[par_height$treatment == "0"] = "control"


par_height$date = as.Date(par_height$date, format="%m.%d.%y")



par_height$PAR_inf[par_height$PAR_inf>21] = round(par_height$PAR_inf[par_height$PAR_inf>21])
par_height$mean_height = rowMeans(par_height[,5:8])
par_height$PAR_diff = par_height$PAR_sup- par_height$PAR_inf

str(par_height)

par_height$species_2[c(par_height$species_2 == "-")] = NA
par_height$species_3[c(par_height$species_3 == "-")] = NA
str(par_height)

ggplot(par_height)+
  facet_wrap(~treatment)+
  geom_boxplot(aes(date, mean_height, col = season))
geom_smooth(aes(date, mean_height), se=FALSE, span = 0.9)
ggplot(par_height)+
  facet_wrap(~treatment)+
  geom_point(aes(date, pheno_sp3, col = season))
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
write.table(par_height, file.path("data","par_and_height.txt"), dec = ".", col.names = TRUE)


#setwd("~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex")

# save(enz_act,file = "organized_data/enzimatic_activity.RData")
# save(par_height,file = "organized_data/par_and_height.RData")
# save(herbi,file = "organized_data/herbivoty_pathogen.RData")
# save(cover,file = "organized_data/cover.RData")
# save(bioma,file = "organized_data/bioma.RData")
# save(interactions, file = "organized_data/interactions_triplet.RData")
nrow(enz_act)
