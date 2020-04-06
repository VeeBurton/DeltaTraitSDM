
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom)
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

?agricolae

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




### use coding club process
# use standardised data
sp3 <- read.csv("./Scots_pine/Scots_pine_H_cent_scal_allvars.csv") # 
sp3$X <- NULL
str(sp3)
colnames(sp3)[1]<-"H"
sp<-na.omit(sp)
sp3$Hraw <- sp$H
sp3$Trial <- sp$Trial
sp3$Provenance <- sp$Provenance

basic.lm <- lm(Hraw ~ Trial:Provenance, data = sp3)
summary(basic.lm)
# plot the data
ggplot(sp3, aes(x =Trial:Provenance, y = H)) +
    geom_jitter() +
    geom_smooth(method = "lm")
# plot residuals
plot(basic.lm, which = 1) # red line should be flat like the dashed grey line
# look at qqplot
plot(basic.lm, which = 2) # points should fall on diagonal line

# check for data independence - essentially checking if need to account for random effects
boxplot(Hraw ~ Trial:Provenance, data = sp3)
# colour plot
ggplot(sp3, aes(x = Trial:Provenance, y = Hraw, colour = Trial:Provenance)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position = "none")
# split plot
ggplot(aes(Trial:Provenance, Hraw), data = sp3) + 
    geom_point() + 
    facet_wrap(~ Trial:Provenance) + # create a facet for each trial
    xlab("Trial:Provenance") + 
    ylab("height")

# include provenance as a fixed effect
prov.lm <- lm(Hraw ~ Trial:Provenance, data = sp3)
tidy(prov.lm)
# as a random effect
mixed.lmer <- lmer(Hraw ~ DD_18_T + (1|Trial/Provenance), data = sp3)
tidy(mixed.lmer)
plot(mixed.lmer) 
qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer)) # points should fall on line


### comparing models (datacamp - intro to statistical modelling)
basic.lm
mixed.lmer
basic_output <- predict(basic.lm, newdata=sp3)
mixed_output <- predict(mixed.lmer, newdata=sp3)
# case-by-case differences
basic_model_differences <- with(sp3, H - basic_output)
mixed_model_differences <- with(sp3, H - mixed_output)
# calculate mean square errors
mean(basic_model_differences ^ 2)
mean(mixed_model_differences ^ 2) # better predictions
