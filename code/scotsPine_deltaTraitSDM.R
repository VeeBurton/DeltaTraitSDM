library(tidyverse)
library(ggplot2)
library(ggmap)
library(lme4)
library(raster)
library(tmap)
library(rgdal)

wd <- "~/R/DeltaTraitSDM/"

# sp data
sp.raw <- read.csv(paste0(wd,"Scots_pine/Scots_pine_H.csv")) # raw data
sp.raw$X <- NULL
str(sp.raw)
sp.raw<-na.omit(sp.raw)

# check data for input to climateEU
cEUinput <- read.csv(paste0(wd,"climate_surfaces/for-climateEU.csv"))
head(cEUinput)
ggplot(cEUinput)+
  geom_point(aes(long,lat))

# use 30 yr average rasters first
climate <- read.csv(paste0(wd,"climate_surfaces/for-climateEU_Normal_1961_1990Y.csv"))
#head(climate)
colnames(climate)<-c('x','y','lat','long','elev','MAT','MWMT','MCMT','TD','MAP','MSP','AHM','SHM','DD0','DD5','DD_18','DD18','NFFD','bFFP','eFFP','FFP','PAS','EMT','Eref','CMD')
climate[climate=="-9999"]<-NA
climate <- na.omit(climate)
summary(climate)
#ggplot(climate)+geom_point(aes(x,y))
coordinates(climate) <- ~ x + y
gridded(climate) <- TRUE
#summary(climate)
MWMT_P <- raster(climate,layer='MWMT')
PAS_T <- raster(climate,layer="PAS")
MAP_T <- raster(climate,layer="MAP")
FFP_P <- raster(climate,layer="FFP")

# process rasters for planting year (T) and when seeds were harvested (P)

# seeds taken (P) | 2007
# planted (T) | 2012

# use Climate EU to generate raster surfaces for provenance and trial sites separately
clim07 <- read.csv(paste0(wd,"climate_surfaces/for-climateEU_Year_2007Y.csv"))
head(clim07)
colnames(clim07)<-c('x','y','lat','long','elev','MAT_P','MWMT_P','MCMT_P','TD_P','MAP_P','MSP_P','AHM_P','SHM_P',
                    'DD0_P','DD5_P','DD_18_P','DD18_P','NFFD_P','bFFP_P','eFFP_P','FFP_P','PAS_P','EMT_P','Eref_P','CMD_P')
clim07[clim07=="-9999"]<-NA
clim07 <- na.omit(clim07)
summary(clim07)

clim12 <- read.csv(paste0(wd,"climate_surfaces/for-climateEU_Year_2012Y.csv"))
head(clim12)
colnames(clim12)<-c('x','y','lat','long','elev','MAT_T','MWMT_T','MCMT_T','TD_T','MAP_T','MSP_T','AHM_T','SHM_T',
                    'DD0_T','DD5_T','DD_18_T','DD18_T','NFFD_T','bFFP_T','eFFP_T','FFP_T','PAS_T','EMT_T','Eref_T','CMD_T')
clim12[clim12=="-9999"]<-NA
clim12 <- na.omit(clim12)
summary(clim12)

coordinates(clim07) <- ~ x + y
gridded(clim07) <- TRUE
MWMT_P <- raster(clim07,layer='MWMT_P')
plot(MWMT_P,main="Mean warmest month temperature at the provenance")
PAS_T <- raster(clim12,layer='PAS_T')
plot(PAS_T,main="Precipitation as snow at the trial site")

# once rasters are processed, scale them, or ensure that model equation scales them
# raster stack of climate variables
ST_H <- MWMT_P
#mean(sp$Hraw) #  stand age
values(ST_H) <- 126 # assign mean height to raster of same extent
plot(ST_H)

rsmax1 <- stack(ST_H, MWMT_P, PAS_T) 
names(rsmax1)= c("Hraw","MWMT_P","PAS_T")
names(rsmax1)
plot(rsmax1)

#crs(predH)<- "+init=epsg:27700" # british national grid
UK <- shapefile(paste0(wd,"Scots_pine/european_region_region.shp"))
unique(UK$NAME)
scot <- UK[which(UK$NAME=="Scotland Euro Region"),]
plot(scot)

pair.mod2 <- lmer(log(W17Height) ~ scale(MWMT_P) + scale(PAS_T) + scale(MWMT_P)*scale(PAS_T) + (1|Trial/Block/Provenance), data=sp.raw)
predH1 <- predict(rsmax1, pair.mod2, re.form=NA, type='response') #whittet.mod2 #dredge1
predH1<- exp(predH1)
predH1_crop <- crop(predH1, scot)
plot(predH1_crop, main="Tree height predicted (pair.mod2)")

#tm_shape(predH1_crop) +
  #tm_style("beaver")+
  #tm_raster(palette="PuBuGn",style = 'pretty') +
  #tm_legend(outside = TRUE)+
  #tm_shape(scot)+
  #tm_borders(lwd=1)+
  #tm_layout(title="Tree height predicted by pair.mod2 (30yr average climate)", title.size = 2)

whittet.mod1 <- lmer(log(W17Height) ~ scale(MAP_T) + scale(FFP_P) + scale(MAP_T)*scale(FFP_P) + (1|Trial/Block/Provenance), data = sp.raw)
rsmax2 <- stack(ST_H, MAP_T, FFP_P) 
names(rsmax2)= c("Hraw","MAP_T","FFP_P")
names(rsmax2)
plot(rsmax2)
predH2 <- predict(rsmax2, whittet.mod1, re.form=NA, type='response') #whittet.mod2 #dredge1
predH2<- exp(predH2)
predH2_crop <- crop(predH2, scot)
plot(predH2_crop, main="Tree height predicted (whittet.mod1)")

# compare overlap of two model outputs
comp <- mean(predH1_crop,predH2_crop)
plot(comp, main="Mean height between 2 models")
# plot side by side with current distribution
nwss <- shapefile(paste0(wd,"Scots_pine/FCPRODUCT_S_NWSS.shp")) 
#summary(nwss)
nwss<-st_as_sf(nwss)
#unique(nwss$DOM_HABITA)
nativePine<-nwss[which(nwss$DOM_HABITA=="Native pinewood"),]
ext <- extent(nativePine)
r <- raster(ext, res=500)
NPraster <- rasterize(nativePine,r,field=1)
plot(NPraster)
stack1 <- stack(comp,NPraster)
plot(stack1)
