
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

library(raster)
coordinates(climate) <- ~ x + y
gridded(climate) <- TRUE
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

rsmax <- stack(ST_H, MWMT_P, PAS_T) 
names(rsmax)= c("Hraw","MWMT_P","PAS_T")
names(rsmax)
plot(rsmax)

pair.mod2 <- lmer(log(Hraw) ~ scale(MWMT_P) + scale(PAS_T) + scale(MWMT_P)*scale(PAS_T) + (1|Trial/Block/Provenance), data=rsmax)
predH <- predict(rsmax, pair.mod2, re.form=NA, type='response') #whittet.mod2 #dredge1
predH<- exp(predH)
plot(predH, main="Tree height predicted")

crs(predH)<- "+init=epsg:27700" # british national grid

tm_shape(predH) +
  tm_style("beaver")+
  tm_raster(palette="PuBuGn",style = 'pretty', title = "Tree height predicted") +
  tm_legend(outside = TRUE)+
  tm_shape(UK)+
  tm_borders(lwd=1)

#whittet.mod1 <- lmer(log(Hraw) ~ MAP_T + FFP_P + MAP_T*FFP_P + (1|Trial/Block/Provenance), data = training)
#rsmax <- stack(ST_H, MAP, FFP) 
#names(rsmax)= c("Hraw","MAP_T","FFP_P")
#names(rsmax)
#plot(rsmax)