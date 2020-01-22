
rm(list=ls())

#Load libraries
library(raster)
library(rgdal) 
library(rworldmap)
library(sjPlot) 
library(lattice)
library(spatial.tools) #libreria para sincronizar rasters ### No
library(maptools)
library(dismo)

### New Extention

newext <- c(-15, 60, 35, 70)

###Trial
#setwd("/home/homero/Documents/Doctorado/DataFagus/RCP8.7_2070")
setwd("C:/Users/vanessa.burton/OneDrive - Forest Research/Documents/R/DeltaTraitSDM-Code")
list.files()

Bio13_Tif<-"bio13_RCP8.5.tif"
Bio13=raster(Bio13_Tif)

Bio13_T <- crop(Bio13, newext)

Bio13_T_Tst <- scale(Bio13_T, center = TRUE, scale = TRUE)
plot(Bio13_T_Tst)


###Provs

Pet.max<-raster("C:/Users/vanessa.burton/OneDrive - Forest Research/Documents/R/DeltaTraitSDM-Code/pet_max_P.gri")

Pet.max_P <- crop(Pet.max, newext)

Pet.max_P_Pst <- scale(Pet.max_P, center = TRUE, scale = TRUE)
plot(Pet.max_P_Pst)


#######  Age  ##########

### Get values

#Import data
#getwd()
#setwd("/home/homero/Documents/Doctorado/R/Scripts_H/Final/Enviar_Marta")
#list.files()
Fagus <- read.table("Fagus_H_Final.txt",header=T)
names(Fagus)
head(Fagus)
### Age 10

FagusAge10<-Fagus[Fagus$Age=="10",]
Age10 <- Bio13_T
unique(FagusAge10$St_Age)
values(Age10) <- 1.464814


#### The model ####
library(lme4)

#the model: 
#setwd("C:/Users/vanessa.burton/OneDrive - Forest Research/Documents/R/DeltaTraitSDM-Code")
M61 <- St_Pet.max_P ~ St_Bio13_T
M61<- lmer(log(H) ~ St_Age + St_Pet.max_P + St_Bio13_T + I(St_Pet.max_P^2) + I(St_Bio13_T^2) + St_Age*St_Pet.max_P + St_Age*St_Bio13_T +  St_Pet.max_P*St_Bio13_T + (1|Trial/Block/Tree_ID) + (1|ID_ProvCode), Fagus)

#load("M61.rda")

M1 <- M61
#summary(M1)
library(tidyverse)
library(broom)
augment(M1)

##### Load Fagus Distribution
getwd()
setwd('/home/homero/Documents/Doctorado/DataFagus/Resultados/Distribution_Fagus_sylvatica')
list.files()

FagusDistr <- shapefile("Fagus_sylvatica_EUFORGEN.shp")




### Load Europe

library("rworldmap")

worldmap <- getMap(resolution = "high")
europe <- worldmap[which(worldmap$REGION=="Europe"),] 

europe <- worldmap[which(grepl("Europe", worldmap$GEO3) & 
                           as.character(worldmap$NAME) != "Russia"),]               
#plot(europe, xlim = c(-25, 50), ylim = c(35, 71), asp = 1)


#### Age 10 ####
St_Age <- Age10
St_Pet.max_P <- Pet.max_P_Pst 
St_Bio13_T <- Bio13_T_Tst
rsmax <- stack(St_Age, St_Pet.max_P, St_Bio13_T)
names(rsmax)=c("St_Age","St_Pet.max_P","St_Bio13_T")
names(rsmax)

predAgemax <- predict(rsmax, M1, re.form = NA, type='response')
#predAgemax <- predHeight(Age=St_Age, var_Pop=St_Pet.max_P, var_Trial=St_Bio13_T)

predAgemax <- exp(predAgemax)

#plot (predAge0, main = "Prediction Height Age 0")

PredAge_max<- mask(predAgemax, FagusDistr) #### Mask

plot(PredAge_max, xlim = c(-10, 45), ylim = c(35, 65))

plot  (europe, add=T)
  


setwd("/home/homero/Documents/Doctorado/DataFagus/Final_Maps/1Trait_Future/Real")
writeRaster(PredAge_max, "PredHeigt_2070_RCP8.5", format = "GTiff")


#### Diference between 2070 and 20015
### 2070
setwd("/home/homero/Documents/Doctorado/DataFagus/Final_Maps/1Trait_Future/Real")
list.files()

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
