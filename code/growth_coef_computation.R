library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(tidyverse)
 
#######################################################
  ########## coef of growth rate per plot ##########
#######################################################

localDir="~/Fac/Cesure2/Plant_species_coexistence/functions_and_coex/"
setwd(localDir)
data_dir = "organized_data"

load(file.path(data_dir,"par_and_height.RData"))
str(par_height)

### Time past since the first measurement (in day)
t0 = as.numeric(as.POSIXct(par_height$date))[1]
par_height$time = (as.numeric(as.POSIXct(par_height$date)) - t0)/(60*60*24)


par_height$plot = as.numeric(as.character(par_height$plot))


par_height_long = par_height%>%
  pivot_longer(cols= c("heigth_A", "heigth_B", "heigth_C", "heigth_D"),
               names_to = "heights_spot",values_to = "height")

ggplot(subset(par_height_long, par_height_long$treatment == "control" ))+
  geom_point(aes(time, height))+
  geom_point(aes(time, mean_height), col ="red")+
  geom_smooth(aes(time, height),formula = y ~ poly(x, 2), method = "lm",
              se=F, linetype = 2, size = 0.5)+

  facet_wrap(~plot)
ggplot(subset(par_height_long, par_height_long$treatment == "nitrogen" ))+
  geom_point(aes(time, height))+
  geom_point(aes(time, mean_height), col ="red")+
  geom_smooth(aes(time, height),formula = y ~ poly(x, 2), method = "lm",
              se=F, linetype = 2, size = 0.5)+
  
  facet_wrap(~plot)

### x²:
par_height$time2 = par_height$time**2

### formula : y = x + x²
formula =  mean_height ~ time + time2

### list of all the models per plot (if needed)
canopy_list_models = list()

### maxtix of all the coefficients per plot
coef_df = matrix(ncol = 5, nrow = max(par_height$plot)*2)

for (n_plot in 1:max(par_height$plot)){
  
  tabl_control = par_height[par_height$plot == n_plot & par_height$treatment == "control",]
  tabl_nitro = par_height[par_height$plot == n_plot & par_height$treatment == "nitrogen",]
  
  
  model_control = lm(formula = formula , data = tabl_control)
  model_nitro = lm(formula = formula , data = tabl_nitro)
  
  canopy_list_models <- append(canopy_list_models,
                               lapply(n_plot, function(x) list( model_control, model_nitro) ))
  
  coef_df[2*n_plot-1,] = c(model_control$coefficients, 0, n_plot )
  coef_df[2*n_plot,] = c(model_nitro$coefficients, 1,n_plot )
}

#names(coef_df) = c("intercept", "coef1", "coef2", "treatment", "plot")

### to match the row of "par_height", for after calculate predictions at each 
# measurement time
temp=coef_df
while(nrow(coef_df) < nrow(par_height)){
  coef_df = rbind(coef_df, temp)
}
nrow(coef_df)
head(coef_df)
predictor = matrix(c(rep(1,nrow(coef_df)),
                     par_height$time,
                     par_height$time2),
                   ncol =3, dimnames = list(c(),c("intercept", "x", "x2")
                   ))
dim(t(predictor))
dim(as.matrix(coef_df[,1:3]))
pred=c()
i=1
### pred = [1,time,time2] %*% t([intecept, coef1, coef2])
# (time2 == time²)
predict(canopy_list_models[[1]][[1]])
i=1
for ( i in 1:nrow(predictor)){
  pred = append(pred, predictor[i,] %*% as.matrix(coef_df[i,1:3]))

}                    
par_height$pred = pred

par_height_long = par_height%>%
  pivot_longer(cols= c("heigth_A", "heigth_B", "heigth_C", "heigth_D"),
               names_to = "heights_spot",values_to = "height")

### check that the calculate prediction match the ggplot prediction
# ( just to see if there is no error in the computation)
ggplot(subset(par_height_long, par_height_long$treatment == "nitrogen" ))+
  geom_point(aes(time, height))+
  geom_point(aes(time, mean_height), col ="red")+
  geom_line(aes(time, pred), col ="red")+
  geom_smooth(aes(time, height),formula = y ~ poly(x, 2), method = "lm",
              se=F, linetype = 2, size = 0.5)+
  
  facet_wrap(~plot)
summary(canopy_list_models[[1]][[2]])
predict(canopy_list_models[[1]][[2]])
tan.poly2 <- function(intercept, coef1, coef2, x1){ 
  #y = f(a) - f'(a)+x
  #y = f(a)+f'(a)(x-a)
  # y = y1 + y'1x - y'1 * x1
  # y = slope + y1-y'1*x1
  y1 = coef2 * x1 ** 2 + coef1 * x1 + intercept
  deriv = 2 * x1 * coef2  + coef1
  
  inter_tang = -x1*deriv + y1
  
  return(c(inter_tang, deriv))
}

###### Compute the tangent line of time = 0 and time = 38 first and last point
tang0 = c()
tang38 = c()
for ( i in 1:nrow(coef_df[1:120,])){
  inter = coef_df[i, 1]
  coef1 = coef_df[i, 2]
  coef2 = coef_df[i, 3]
  tang0 = rbind(tang0,tan.poly2(inter, coef1, coef2, 0))
  tang38 = rbind(tang38,tan.poly2(inter, coef1, coef2, 38))
  
} 

canopy_height_growth = data.frame(tang0_intecept = tang0[,1], tang0_slope = tang0[,2],
           tang38_intecept = tang38[,1], tang38_slope = tang38[,2],
           plot_ID = par_height$plot_ID[1:120],
           treatment = par_height$treatment[1:120])

par_height$tang0_intecept = tang0[,1]
par_height$tang0_slope = tang0[,2]
par_height$tang38_intecept = tang38[,1]
par_height$tang38_slope = tang38[,2]
str(par_height_long)

summary(canopy_list_models[[9]][[1]])
ggplot(subset(par_height_long, par_height_long$treatment == "nitrogen"&
              par_height_long$plot_ID ==34 ))+
  geom_point(aes(time, height), cex =2)+
  geom_line(aes(time, pred), col ="red", size =1)+
  facet_wrap(~plot)+
  geom_abline(data = subset(par_height, par_height$treatment == "nitrogen"&
                            par_height$plot_ID ==34  ),
              aes(intercept = tang0_intecept, slope = tang0_slope),
              col = "green", linetype =2, size =1)+
  geom_vline(aes(xintercept = 38 ), col = "blue", size =1, linetype =2)+
  xlab("time (in days)") + ylab("Canopy height")
  #geom_abline(data = subset(par_height, par_height$treatment == "nitrogen" ),
  #         aes(intercept = tang38_intecept, slope = tang38_slope), col = "blue")
max_canopy_heigth = c()
max_canopy_time = c()
for(n_plot in 1:60){
  max_canopy_heigth = append(max_canopy_heigth,
                             max(par_height[c(par_height$plot_ID == n_plot &
                                            par_height$treatment == "control"),]$pred))
  index = which.max(par_height[c(par_height$plot_ID == n_plot &
                     par_height$treatment == "control"),]$pred) 
  max_canopy_time = append(max_canopy_time,
                           par_height[c(par_height$plot_ID == n_plot &
                            par_height$treatment == "control"),]$time[index])
  
  
  max_canopy_heigth = append(max_canopy_heigth,
                             max(par_height[c(par_height$plot_ID == n_plot &
                                            par_height$treatment == "nitrogen"),]$pred))
  index = which.max(par_height[c(par_height$plot_ID == n_plot &
                                   par_height$treatment == "nitrogen"),]$pred) 
  max_canopy_time = append(max_canopy_time,
                           par_height[c(par_height$plot_ID == n_plot &
                                          par_height$treatment == "nitrogen"),]$time[index])
  
}
max_canopy_time
tail(par_height)
names(par_height_long)
canopy_height_growth$max_canopy_heigth = max_canopy_heigth
canopy_height_growth$max_canopy_time = max_canopy_time


save(par_height,file = "organized_data/par_and_height.RData")  
save(canopy_height_growth, file = "organized_data/canopy_height_growth.RData")
