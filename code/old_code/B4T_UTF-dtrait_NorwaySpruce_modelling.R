###---Trial data
#-reference climate for trial data
library(raster)
clmTR<-stack(list.files('/media/maurizio/mauHD/GIS_generali/ClimateEU_generated/Normal_81-10_1km_Albers/',full.names=T,pattern='.tif'))
names(clmTR)
plot(clmTR[['MAT']])
trData<-read.csv('/media/maurizio/mauHD/lavoro/CNR/B4EST/WP1_ZhiqiangChen_data/IUFRO196468_trials_adj.csv')
head(trData)
trData<-SpatialPointsDataFrame(coords=trData[,c('Lng','Lat')],proj4string=CRS('+init=epsg:4326'),data=trData)
trData<-spTransform(trData,crs(clmTR))
points(trData,pch=19,cex=0.5)
trData<-extract(clmTR,trData,sp=T)
summary(trData)
trData<-trData@data[-which(is.na(trData@data[,'AHM'])),]
head(trData)
trDataFULL<-NULL
for(i in 1:nrow(trData)) {
  #i<-1
  trDataFULL<-rbind(trDataFULL,cbind(data.frame(Site=paste(trData[i,'Site'],'_',which(trData[i,7:17]=='x'),sep='')),trData[i,names(clmTR)]))
}
nrow(trDataFULL)
head(trDataFULL)

###---Norway Spruce data
#-reference climate
clmNS<-stack(list.files('/media/maurizio/mauHD/GIS_generali/ClimateEU_generated/Decadal_1961-1970_forB4EST/',full.names=T,pattern='.tif'))
plot(clmNS[['MAT']])
nsData<-read.csv('/media/maurizio/mauHD/lavoro/CNR/B4EST/WP1_ZhiqiangChen_data/IUFRO196468_provenances_adj.csv')
nsData<-subset(nsData,LAT>0 & LNG>0)
head(nsData)
summary(nsData)
nsDataSP<-SpatialPointsDataFrame(coords=nsData[,c('LNG','LAT')],proj4string=CRS('+init=epsg:4326'),data=nsData)
nsDataSP<-spTransform(nsDataSP,crs(clmNS))
points(nsDataSP,pch=19,cex=0.5)

#-sdm to mask the mxm prediction
# library(dismo)
# mn<-maxent(clmNS,nsDataSP,nbg=1100)
# #mn
# plot(mn)
# sdm<-predict(mn,clmTR)
# writeRaster(sdm,'/media/maurizio/mauHD/lavoro/CNR/B4EST/WP1_ZhiqiangChen_data/NSsdm.tif',format='GTiff',overwrite=T)
sdm<-raster('/media/maurizio/mauHD/lavoro/CNR/B4EST/WP1_ZhiqiangChen_data/NSsdm.tif')
plot(sdm,main='SDM')
points(nsDataSP,pch=19,cex=0.5)

nsData<-extract(clmNS,nsDataSP,sp=T)
summary(nsData)
nsData<-nsData@data[-which(is.na(nsData@data[,'MAT'])),]
head(nsData)
nsDataTEMP<-NULL
for(i in 1:nrow(nsData)) {
  #i<-1
  a<-nsData[i,c(3:4,9:ncol(nsData))]
  nsDataTEMP<-rbind(nsDataTEMP,cbind(data.frame(ProvID=a$ProvID,Site=paste(names(a)[3:16],a[,'Block'],sep='_'),HT=as.numeric(a[,3:16])),a[,17:ncol(a)]))
}
nrow(nsDataTEMP)
head(nsDataTEMP)

###---let's see befor merging where Trials were established
plot(MAP~MAT,nsDataTEMP,pch=19,xlab='Mean annual Temperature (°C)',ylab='Mean toatl annual precipitation (mm)',main='All test sites are inside the range')
points(MAP~MAT,trDataFULL,pch=19,col='red')
spPCA<-nsDataTEMP[,names(clmNS)]
row.names(spPCA)<-paste('Pr',1:nrow(spPCA),sep='_')
trPCA<-trDataFULL[,names(clmNS)]
row.names(trPCA)<-paste('Tr',1:nrow(trPCA),sep='_')
dbPCA<-rbind(spPCA,trPCA)
pca<-prcomp(dbPCA,center=T,scale=T)
summary(pca)
biplot(pca)

###---mixed modelling analysis
nsDataALL<-merge(nsDataTEMP,trDataFULL[c('Site',names(clmTR))],by='Site')
summary(nsDataALL)
nsDataALL<-nsDataALL[-which(is.na(nsDataALL$HT)),]
names(nsDataALL)
summary(nsDataALL)
###--looking for relationships between HT and the provenance climate minus site climate
# y<-nsDataALL[,'HT']
# predictionDF<-data.frame(x=seq(-100000,10000,1))
# for(i in 4:ncol(nsDataALL)) {
#   #i<-5
#   x<-nsDataALL[,names(nsDataALL)[i]]
#   plot(y~x,ylab='Standardised height',xlab=names(nsDataALL)[i],pch=19)
#   library(mgcv)
#   LM<-lm(y~x)
#   GM<-gam(y~s(x))
#   lines(predict(LM,predictionDF)~predictionDF[,1],col='red',lwd=2)
#   lines(predict(GM,predictionDF)~predictionDF[,1],col='blue',lwd=2)
# }
###---data scaling
nsDataALLscaled<-cbind(nsDataALL[,1:3],scale(nsDataALL[,4:ncol(nsDataALL)]))
library(agricolae)
a<-correlation(cbind(nsDataALLscaled$HT,nsDataALLscaled[,4:23]),method='spearman') #climate of the Provenance
#a<-correlation(cbind(nsDataALLscaled$HT,nsDataALLscaled[,24:43]),method='spearman') #climate of the Site
#a<-correlation(cbind(nsDataALLscaled$HT,nsDataALLscaled[,24:43]-nsDataALLscaled[,4:23]),method='spearman') #difference
a$correlation
library(lmerTest)
lmx<-lmer(log(HT)~AHM.x+AHM.y+AHM.x*AHM.y+PAS.x+PAS.y+PAS.x*PAS.y+FFP.x+FFP.y+FFP.x*FFP.y+(1|Site),data=nsDataALLscaled,REML=T)
summary(lmx)
plot(lmx,pch=19,cex=0.5,xlab='Fitted values',ylab='Residuals')
library(MuMIn)
r.squaredGLMM(lmx,pj2014=T)
#-prediction
rsRAST<-scale(stack(clmNS[['AHM']],clmTR[['AHM']],clmNS[['PAS']],clmTR[['PAS']],clmNS[['FFP']],clmTR[['FFP']]))
names(rsRAST)<-names(lmx@frame)[2:I(ncol(lmx@frame)-1)]
VmapMXMprd<-exp(predict(rsRAST,lmx,re.form=NA,type='response'))
plot(VmapMXMprd)
Msdm<-sdm>0.47
values(Msdm)[which(values(Msdm)==0)]<-NA
plot(Msdm)
plot(sdm,col='grey95',legend=F)
plot(mask(VmapMXMprd,Msdm),add=T)

###---Generalised Additive Mixed-effects Model 
library(gamm4)
gmx<-gamm4(HT~s(AHM.x)+s(AHM.x,by=AHM.y)+s(PAS.x)+s(PAS.x,by=PAS.y)+s(FFP.x)+s(FFP.x,by=FFP.y),random=~(1|Site),
           data=nsDataALLscaled,REML=T,verbose=T) #[sample(1:nrow(nsDataALLscaled),size=250),]
summary(gmx$mer)
summary(gmx$gam)
par(mfrow=c(3,2))
plot(gmx$gam)
par(mfrow=c(1,1))
VmapMXMprd<-predict(rsRAST,gmx$gam,re.form=NA,type='response')
#pdf('/media/maurizio/mauHD/lavoro/CNR/B4EST/WP1_ZhiqiangChen_data/NorwaySpruce_GAMM_prediction.pdf',height=20,width=20)
plot(VmapMXMprd,zlim=c(0,120),main='Spatial prediction for Norway spruce data using Generalised Additive Mixed Models (GAMM)')
Msdm<-sdm>0.47
values(Msdm)[which(values(Msdm)==0)]<-NA
plot(Msdm)
plot(sdm,col='grey95',legend=F)
plot(mask(VmapMXMprd,Msdm),add=T)
dev.off()
