localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "organized_data"

load(file.path(data_dir,"bioma.RData"))
load(file.path(data_dir,"coexistance.RData"))
load(file.path(data_dir,"cover.RData"))
load(file.path(data_dir,"enzimatic_activity.RData"))
load(file.path(data_dir,"herbivoty_pathogen.RData"))
load(file.path(data_dir,"par_and_height.RData"))

str(bioma)
bioma$species = apply( bioma[,c("species_1", "species_2", "species_3")], 1, 
                       function(x) paste(x[!is.na(x)], collapse = "_"))  

bioma$nitrogen = as.character(bioma$nitrogen)
bioma$nitrogen[bioma$nitrogen == "1"] = "nitrogen"
bioma$nitrogen[bioma$nitrogen == "0"] = "control"
names(bioma)[names(bioma) == "plot" | names(bioma) == "nitrogen"] = c("treatment", "plot_ID")

str(cover)

cover$nitrogen[cover$nitrogen == "1"] = "nitrogen"
cover$nitrogen[cover$nitrogen == "0"] = "control"
names(cover)[names(cover) == "plot" | names(cover) == "nitrogen"] = c("plot_ID", "treatment")

str(enz_act)

enz_act$nitrogen = as.character(enz_act$nitrogen)
enz_act$nitrogen[enz_act$nitrogen == "1"] = "nitrogen"
enz_act$nitrogen[enz_act$nitrogen == "0"] = "control"
names(enz_act)[names(enz_act) == "plot" | names(enz_act) == "nitrogen"] = c("treatment", "plot_ID")

str(herbi)

herbi$treatement = as.character(herbi$treatement)
herbi$treatement[herbi$treatement == "1"] = "nitrogen"
herbi$treatement[herbi$treatement == "0"] = "control"
names(herbi)[names(herbi) == "plot" ] = "plot_ID"
herbi$plot_ID = as.character(herbi$plot_ID)

str(par_height)

par_height$nitrogen = as.character(par_height$nitrogen)
par_height$nitrogen[par_height$nitrogen == "1"] = "nitrogen"
par_height$nitrogen[par_height$nitrogen == "0"] = "control"
names(par_height)[names(par_height) == "plot" | names(par_height) == "nitrogen"] = c("plot_ID", "treatment")

save(enz_act,file = "organized_data/enzimatic_activity.RData")
save(par_height,file = "organized_data/par_and_height.RData")
save(herbi,file = "organized_data/herbivoty_pathogen.RData")
save(cover,file = "organized_data/cover.RData")
save(bioma,file = "organized_data/bioma.RData")
