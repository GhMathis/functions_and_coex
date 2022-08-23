localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "oster_data_2022"

##### code to merge control coexistence data and nitrogen coexistence data in one df
coexistance_N <- read.delim2(file.path(data_dir,"biomcoexistence_all_2020_nitrogen_18sp2-3.txt"),
             header =TRUE, dec =".")

coexistance_C <- read.delim2(file.path(data_dir,"final_triple_selection_corrected.txt"),
                            header =TRUE, dec =".")
str(coexistance_N)
str(coexistance_C)

coexistance = coexistance_N[coexistance_N$species %in% coexistance_C$species,]
names(coexistance_N)
names(coexistance_C)
coexistance_C = coexistance_C[,1:10]
coexistance = rbind(coexistance, coexistance_C)
head(coexistance)
str(coexistance)

save(coexistance, file = "organized_data/coexistance.RData")
