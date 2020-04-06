
library(tidyverse)
library(lme4)
library(lmerTest)
library(agricolae)
library(viridis)

#################

sp <- read.csv("./Scots_pine/Scots_pine_H.csv") # raw data
sp$X <- NULL
str(sp)
colnames(sp)[11]<-"H"

#boxplot(H~Trial:Provenance,sp)
sp %>% 
  ggplot(aes(Trial:Provenance,H, colour=Trial))+geom_boxplot()+
  theme_minimal()+
  ylab("Height (mm)")+
  theme(axis.text.x = element_text(size=6,angle = 90, hjust = 1))

#boxplot(H~Provenance,sp)
sp %>% 
  ggplot(aes(Provenance,H, colour=Provenance))+geom_boxplot()+
  theme_minimal()+
  ylab("Height (mm)")+
  theme(axis.text.x = element_blank())

y<-sp$H
x<-paste(sp$Provenance,sp$Trial,sep='_')
m<-lm(y~x) # model of height by trial:provenance
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

# mean height per provenance within each trial site
mean_H<-aggregate(H~Trial*Provenance,sp2,FUN=mean)
#sp2<-na.omit(sp2)
sp2<- sp2 %>% 
  group_by(Trial,Provenance) %>% 
  mutate(mean_H = mean(H))
sp2$mean_H

# all variables from PC1 of Richard's thesis
# growing degree-days, monthly mean temps for Feb and July, annual precipitation, extreme temperature range

# growing degree days
sp2<- sp2 %>% 
  group_by(Trial,Provenance) %>% 
  mutate(mean_DD18P = mean(DD_18_P),
         mean_DD18T = mean(DD_18_T),
         mean_DD18d = mean(DD18diff))

sp2 %>% 
  ggplot(aes(mean_DD18T,H,colour=Provenance))+geom_point()+stat_smooth()

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
