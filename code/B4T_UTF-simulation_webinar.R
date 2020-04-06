####----UTF simulation  ---->>>> each simulation will generate a different output
wd<-'/home/maurizio/Documents/UK_climate/'
#wd<-'/media/maurizio/mauHD/GIS_generali/ClimateEU_generated/UK_only/'

##---loading raster map for spatialization
library(raster)
mcmt<-raster(paste(wd,'Baseline/MCMT.tif',sep=''))
plot(mcmt)

#-suppose to have *n* locations to be used as both provenance and trial site
n<-100
pts<-sampleRandom(mcmt,size=n,sp=T)
points(pts,pch=19,cex=0.7)

#-and that all the n provenances are able to achieve a value of 350 m3.ha-1
vols<-data.frame(ID=1:n,Vol=rep(350,n),Coef=NA,Vol2=NA)
#head(vols)

# then this value will be artificially reduced according to the "distance" between the MCMT of the site and the MCMT of the provenance(s) to generate a matrix of 100*100 values
# with each row as a provenance (average value) and each column as test site
volData<-data.frame(Provenance=paste('P',1:n,sep=''))
for(i in 1:nrow(pts@data)) {
  #i<-1
  #mlt<-sort(sample(seq(0.4,1.1,0.02),n/2,replace=T))
  mlt<-sort(rnorm(n,mean=1,sd=0.4),decreasing=T)
  vols[order(abs(pts@data[,'MCMT']-pts@data[i,'MCMT'])),'Coef']<-mlt
  vols[,'Vol2']<-vols[,'Vol']*vols[,'Coef']
  volData<-cbind(volData,data.frame(Volumes=vols[,'Vol2']*rnorm(n,mean=1,sd=0.2)))
  #plot(volData[,i+1]~I(pts@data[,'MCMT']-pts@data[i,'MCMT']),pch=19)
}
names(volData)<-c('Provenance',paste('Site',1:n,sep='_'))
head(volData)

###---A B and C values for each test site
x<-pts@data[,'MCMT']
ABCval<-NULL
for(i in 2:ncol(volData)) {
  y<-volData[,i]
  fnls<-nls(y~A/(1+((x-B)/C)^2),start=list(A=300,B=3,C=2))
  summary(fnls)
  # #visual check of the fitting
  # 1-deviance(fnls)/deviance(lm(y~1))
  # plot(y~x,pch=19)
  # lines(predict(fnls,data.frame(x=seq(-10,10,0.01)))~seq(-10,10,0.01),col='red',lwd=2)
  ABCval<-rbind(ABCval,data.frame(Site=i-1,A=coef(fnls)['A'],B=coef(fnls)['B'],C=coef(fnls)['C'],R2=1-deviance(fnls)/deviance(lm(y~1))))
}
summary(ABCval)

###---A-B-C raster surfaces (this is the weakest part because it is difficult to generate an artificial dataset to be easily modelled automatically)
#the correlation matrix can be used to derite the most suitable climatic parameter
clm<-stack(list.files(paste(wd,'Baseline',sep=''),full.names=T,pattern='.tif'))
aa<-extract(clm,pts,df=T)
library(agricolae)
crl<-correlation(cbind(ABCval[,2:4],aa[,-1]),method='spearman')
ABCval<-cbind(ABCval,aa[,-1])

##---we will use GAM because easier and much more flexible than LM
library(mgcv)
sort(crl$correlation[1,])
Afn<-gam(A~s(AHM)+s(PPT_wt),data=ABCval)
summary(Afn)
sort(crl$correlation[2,])
Bfn<-gam(B~s(MCMT),data=ABCval)
summary(Bfn)
sort(crl$correlation[3,])
Cfn<-gam(C~s(EMNT)+s(DD_0),data=ABCval)
summary(Cfn)

#-caculating raster maps (quite slow don't worry...we will probably reduce the resolution)
#clm<-aggregate(clm,fact=4,fun=mean) #from 250m to 1km
Amap<-predict(clm,Afn)
Bmap<-predict(clm,Bfn)
Cmap<-predict(clm,Cfn)
#-Final map of predicted Volumes for current time (1961-1990)
Vmap<-Amap/(1+((mcmt-Bmap)/Cmap)^2)
plot(Vmap)
#-Final map of predicted Volumes for a future scenaros (RCP4.5 and RCP8.5 2050s)
Fmcmt<-stack(list.files(paste(wd,'/Future/',sep=''),full.names=T))
VmapFUTURE<-Amap/(1+((Fmcmt-Bmap)/Cmap)^2)
names(VmapFUTURE)<-c('RCP4.5_2050s','RCP8.5_2050s')
plot(VmapFUTURE) #prediction
plot(Vmap-VmapFUTURE,main=c('RCP4.5_2050s','RCP8.5_2050s')) #-variation

##
####
##########----what about with a blackbox algorithm such as RF? Much higher expl.Ver but much slower
ABCval<-cbind(ABCval[,-6],aa[,-1])
head(ABCval)
library(randomForest)
Afn<-randomForest(A~.,data=ABCval[,c(2,6:ncol(ABCval))],ntree=5000)
Afn
Bfn<-randomForest(B~.,data=ABCval[,c(3,6:ncol(ABCval))],ntree=5000)
Bfn
Cfn<-randomForest(C~.,data=ABCval[,c(4,6:ncol(ABCval))],ntree=5000)
Cfn
AmapRF<-predict(clm,Afn)
BmapRF<-predict(clm,Bfn)
CmapRF<-predict(clm,Cfn)
#-Final map of predicted Volumes for current time (1961-1990) with RF A-B-C maps
VmapRF<-AmapRF/(1+((mcmt-BmapRF)/CmapRF)^2)
plot(VmapRF)
#-Final map of predicted Volumes for a future scenaros (RCP4.5 and RCP8.5 2050s)
Fmcmt<-stack(list.files(paste(wd,'/Future/',sep=''),full.names=T))
VmapFUTURErf<-AmapRF/(1+((Fmcmt-BmapRF)/CmapRF)^2)
names(VmapFUTURErf)<-c('RCP4.5_2050s','RCP8.5_2050s')
plot(VmapFUTURErf) #prediction
plot(VmapRF-VmapFUTURErf,main=c('RCP4.5_2050s','RCP8.5_2050s')) #-variation

##
####
##########----DeltaTraitSDM method
head(volData)
data4mxm<-stack(volData)
names(data4mxm)<-c('Volume','ProvTrial')
data4mxm$Provenance<-paste('P',1:100,sep='')
data4mxm$MCMT_provTr<-rep(pts@data$MCMT,each=n)
data4mxm$MCMT_origin<-rep(pts@data$MCMT,n)
head(data4mxm)
library(lmerTest)
mxm<-lmer(Volume~MCMT_provTr+MCMT_origin+MCMT_provTr*MCMT_origin+(1|ProvTrial),data=data4mxm,REML=T)
summary(mxm)
plot(mxm)
1-deviance(mxm)/deviance(lm(data4mxm$Volume ~1))
VmapMXM<-mean(coef(mxm)$ProvTrial[,1]+mcmt*(mean(coef(mxm)$ProvTrial[,2])+mean(coef(mxm)$ProvTrial[,3]+mean(coef(mxm)$ProvTrial[,4]))))
plot(VmapMXM)
###---end