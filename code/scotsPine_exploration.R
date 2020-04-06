
library(tidyverse)
library(lme4)
library(lmerTest)
library(agricolae)

#################

sp <- read.csv("./Scots_pine/Scots_pine_H.csv") # raw data
sp$X <- NULL
str(sp)

boxplot(W17Height~Trial:Provenance,sp)
boxplot(W17Height~Provenance,sp)
y<-sp$W17Height
x<-paste(sp$Provenance,sp$Trial,sep='_')
m<-lm(y~x)
summary(m)
anova(m)
a<-duncan.test(m,trt="x")
plot(a) # this suggests there are significant differences between provenance at trial sites

##################

sp2 <- sp %>% mutate(MATdiff = MAT_P-MAT_T,
                     MWMTdiff = MWMT_P-MWMT_T,
                     MCMTdiff = MCMT_P-MCMT_T,
                     TDdiff = TD_P-TD_T,
                     MAPdiff = MAP_P-MAP_T,
                     MSPdiff = MSP_P-MSP_T,
                     AHMdiff = AHM_P-AHM_T,
                     SHMdiff = SHM_P-SHM_T,
                     DD0diff = DD0_P-DD0_T,
                     DD18diff = DD_18_P-DD_18_T,
                     NFFDdiff = NFFD_P-NFFD_T,
                     bFFPdiff = bFFP_P-bFFP_T,
                     eFFPdiff = eFFP_P-eFFP_T,
                     FFPdiff = FFP_P-FFP_T,
                     PASdiff = PAS_P-PAS_T,
                     EMTdiff = EMNT_P-EMNT_T,
                     Erefdiff = Eref_P-Eref_T,
                     CMDdiff = CMD_P-CMD_T)

head(sp2[,c(56:73)])
colnames(sp2)[11]<-"H"

# mean height per provenance within each trial site
mean_H<-aggregate(H~Trial*Provenance,sp2,FUN=mean)

# all variables from PC1 of Richard's thesis
# growing degree-days, monthly mean temps for Feb and July, annual precipitation, extreme temperature range

# growing degree days
DD18diff<-aggregate(DD18diff~Trial*Provenance,sp2,FUN=mean)
plot(mean_H[,3]~DD18diff[,3],col=mean_H$Provenance,pch=19) # no pattern/relationship
DD18p<-aggregate(DD_18_P~Trial*Provenance,sp2,FUN=mean)
plot(mean_H[,3]~DD18p[,3],col=mean_H$Provenance,pch=19) # no pattern/relationship
DD18t<-aggregate(DD_18_T~Trial*Provenance,sp2,FUN=mean)
plot(mean_H[,3]~DD18t[,3],col=mean_H$Provenance,pch=19) # inverse bell curve?

# annual precipitation
PASdiff<-aggregate(PASdiff~Trial*Provenance,sp2,FUN=mean)
plot(mean_H[,3]~PASdiff[,3],col=mean_H$Provenance,pch=19)
PASp<-aggregate(PAS_P~Trial*Provenance,sp2,FUN=mean)
plot(mean_H[,3]~PASp[,3],col=mean_H$Provenance,pch=19) # no pattern/relationship
PASt<-aggregate(PAS_T~Trial*Provenance,sp2,FUN=mean)
plot(mean_H[,3]~PASt[,3],col=mean_H$Provenance,pch=19) # inverse bell curve?
