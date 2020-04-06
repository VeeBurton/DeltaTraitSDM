
#rm(list=ls())

#Load libraries
library(raster)
library(rgdal) 
library(rworldmap)
library(sjPlot) 
library(lattice)
library(spatial.tools) 
library(maptools)
library(dismo)

### new extent
newext <- c(-15, 60, 35, 70)

### Trial
setwd("C:/Users/vanessa.burton/OneDrive - Forest Research/Documents/R/DeltaTraitSDM")
list.files()

# so, max temp of warmest month found to be most important explanatory variable for the trial sites?
# or is it max precipitation of wettest month (biovars bio13) (yes, confirmed with marta)

# current climate
Bio13_cc<-raster('./Fagus/bio13.gri')
Bio13_cc<-crop(Bio13_cc,newext)
Bio13_cc <- scale(Bio13_cc, center = TRUE, scale = TRUE) # scale - centre subtracts layer means, scale divides layers by standard deviations
plot(Bio13_cc,  main="Bio 13 (current climate)")

# future climate
Bio13_8.5<-"./Fagus/bio13_RCP8.5.tif" 
Bio13_8.5<-raster(Bio13_8.5)
Bio13_8.5 <- crop(Bio13_8.5, newext)
Bio13_8.5 <- scale(Bio13_8.5, center = TRUE, scale = TRUE) 
plot(Bio13_8.5, main='Bio 13 (RCP8.5)')


### Provenance

# max precipitation most important variable for the provenance sites?

Pet.max<-raster("./Fagus/pet_max_P.gri")
Pet.max<- crop(Pet.max, newext)
Pet.max<- scale(Pet.max, center = TRUE, scale = TRUE)
plot(Pet.max, main='Maximum precipitation (current climate)')


#######  Trait data - Age  ##########

library(tidyverse)

### Get values

#Import data
Fagus <- read.table("./Fagus/Fagus_H_Final.txt",header=T)
#names(Fagus)
#head(Fagus)
str(Fagus)

### Age 10
FagusAge10<-filter(Fagus, Age==10)
Age10 <- Bio13_cc
unique(FagusAge10$St_Age) #  stand age?
values(Age10) <- 1.464814 # assign stand age to raster of bio13?


#### The model ####
library(lme4)

#the model: 
#setwd("C:/Users/vanessa.burton/OneDrive - Forest Research/Documents/R/DeltaTraitSDM-Code")
M61 <- St_Pet.max_P ~ St_Bio13_T
M61<- lmer(log(H) ~ St_Age + St_Pet.max_P + St_Bio13_T + I(St_Pet.max_P^2) + I(St_Bio13_T^2) + St_Age*St_Pet.max_P + St_Age*St_Bio13_T +  St_Pet.max_P*St_Bio13_T + (1|Trial/Block/Tree_ID) + (1|ID_ProvCode), Fagus)

#load("M61.rda")

M1 <- M61
summary(M1)
#library(tidyverse)
#library(broom)
#augment(M1)

##### Load Fagus Distribution
#getwd()
#setwd('/home/homero/Documents/Doctorado/DataFagus/Resultados/Distribution_Fagus_sylvatica')
#list.files()

FagusDistr <- shapefile("./Fagus/Fagus_sylvatica_distribution/Fagus_sylvatica_EUFORGEN.shp") 

### Load Europe

library("rworldmap")

worldmap <- getMap(resolution = "high")
europe <- worldmap[which(worldmap$REGION=="Europe"),] 
europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & as.character(worldmap$NAME) != "Russia"),]               
plot(europe, xlim = c(-25, 50), ylim = c(35, 71), asp = 1)

#### Age 10 ####
St_Age <- Age10
St_Pet.max_P <- Pet.max 
St_Bio13_T <- Bio13_cc
rsmax <- stack(St_Age, St_Pet.max_P, St_Bio13_T) # stack of stand age, max precip, and max temp of warmest month
names(rsmax)=c("St_Age","St_Pet.max_P","St_Bio13_T")
names(rsmax)

predAge10 <- predict(rsmax, M1, re.form=NA, type='response')
predH_Age10<- exp(predAge10)
plot(predAge10, main="Tree height predicted at age = 10")

#### Mask
PredHAge_10<- mask(predH_Age10, FagusDistr) 
plot(PredHAge_10, xlim = c(-10, 45), ylim = c(35, 65), main="Tree height prediction at age = 10 within beech range for current climatic conditions")
plot(europe, add=T)

# and for future climate
St_Age <- Age10
St_Pet.max_P <- Pet.max 
St_Bio13_T <- Bio13_8.5
rsmax <- stack(St_Age, St_Pet.max_P, St_Bio13_T) # stack of stand age, max precip, and max temp of warmest month
names(rsmax)=c("St_Age","St_Pet.max_P","St_Bio13_T")
names(rsmax)

predAge10 <- predict(rsmax, M1, re.form=NA, type='response')
predH_Age10<- exp(predAge10)
plot(predAge10, main="Tree height predicted at age = 10")

#### Mask
PredHAge_10<- mask(predH_Age10, FagusDistr) 
plot(PredHAge_10, xlim = c(-10, 45), ylim = c(35, 65), main="Tree height prediction at age = 10 within beech range for future climatic conditions")
plot(europe, add=T)

#predAgemax <- predict(rsmax, M1, re.form = NA, type='response') # predict height based on future precipitation data
#predAgemax <- predHeight(Age=St_Age, var_Pop=St_Pet.max_P, var_Trial=St_Bio13_T)
#predAgemax <- exp(predAgemax)
#plot (predAge0, main = "Prediction Height Age 0")
#PredAge_max<- mask(predAgemax, FagusDistr) #### Mask
#plot(PredAge_max, xlim = c(-10, 45), ylim = c(35, 65))
#plot(europe, add=T)
  
#setwd("/home/homero/Documents/Doctorado/DataFagus/Final_Maps/1Trait_Future/Real")
#writeRaster(PredAge_max, "PredHeigt_2070_RCP8.5", format = "GTiff")


#### Diference between 2015 and 2070
### 2070
#setwd("/home/homero/Documents/Doctorado/DataFagus/Final_Maps/1Trait_Future/Real")
#list.files()

Height_2070_Tif <- "PredHeigt_2070_RCP8.tif" 
Height_2070 <- raster(Height_2070_Tif)
plot(Height_2070)

### 2015
setwd("/home/homero/Documents/Doctorado/DataFagus/Final_Maps/1Trait_Present")
list.files()
Height_2015_Tif <- "PredHeight_Age10.tif" 
Height_2015 <- raster(Height_2015_Tif)
plot(Height_2015)

#### Difference
Height_Diff <- Height_2015
values(Height_Diff) <- values(Height_2070) - values(Height_2015)
plot(Height_Diff)
setwd("/home/homero/Documents/Doctorado/DataFagus/Final_Maps/1Trait_Future/Diff")
writeRaster(Height_Diff, "PredHeigt_Diff", format = "GTiff")
list.files()

Height_Diff <- "PredHeigt_Diff.tif" 
Height_Diff <- raster(Height_Diff)
plot(Height_Diff)
hist(Height_Diff)
quantile(Height_Diff, prob = c(0.01, 0.99))
H <- Height_Diff
plot(H)
hist(H)

values(H)[values(H) > 119] <- 119
values(H)[values(H) < -33] <- -33
plot(H)
setwd("/home/homero/Documents/Doctorado/DataFagus/Final_Maps/1Trait_Future/Diff")
writeRaster(H, "PredHeigt_Diff_Final", format = "GTiff")
