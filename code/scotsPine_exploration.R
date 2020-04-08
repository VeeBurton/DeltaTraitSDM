
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
a<-agricolae::duncan.test(m,trt="x")
plot(a) # this suggests there are significant differences between provenance at trial sites
head(a)

# try with other variables
y<-sp$bFFP_P
x<-paste(sp$Provenance,sp$Trial,sep='_')
m<-lm(y~x) # model of height by trial:provenance
summary(m)
anova(m)
b<-agricolae::duncan.test(m,trt="x")
plot(b) # this suggests there are significant differences between provenance at trial sites
head(b)

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
#mean_H<-aggregate(H~Trial:Provenance,sp2,FUN=mean)
#sp2<-na.omit(sp2)

# all variables from PC1 of Richard's thesis
# growing degree-days, monthly mean temps for Feb and July, annual precipitation, extreme temperature range
mean_diffs <- sp2 %>% 
  group_by(Trial,Provenance) %>% 
  summarise(meanH = mean(H, na.rm=TRUE),
            mean_DD18p = mean(DD_18_P),
            mean_DD18t = mean(DD_18_T),
            mean_DD18d = mean(DD18diff),
            mean_PASp = mean(PAS_P),
            mean_PASt = mean(PAS_T),
            mean_PASd = mean(PASdiff),
            mean_EMTp = mean(EMNT_P),
            mean_EMTt = mean(EMNT_T),
            mean_EMTd = mean(EMTdiff),
            mean_MCMTp = mean(MCMT_P),
            mean_MCMTt = mean(MCMT_T),
            mean_MCMTd = mean(MCMTdiff),
            mean_MWMTp = mean(MWMT_P),
            mean_MWMTt = mean(MWMT_T),
            mean_MWMTd = mean(MWMTdiff),
            mean_TDp = mean(TD_P),
            mean_TDt = mean(TD_T),
            mean_TDd = mean(TDdiff))

mean_diffs %>% 
  gather(key='variable',value='mean', -Trial,-Provenance,-meanH) %>% 
  filter(variable %in% c("mean_DD18p","mean_DD18t","mean_EMTd","mean_MCMTd","mean_MWMTd","mean_TDd")) %>% 
  ggplot(aes(mean,meanH))+
    facet_wrap(~variable)+
    geom_point()+
    stat_smooth()

# look at variables with largest differences between provenances and trials
summary(sp2[,c(56:73)])
# MAP, MSP, DD0, DD18, NFFD, bFFP, eFFP, FFP, PAS, Eref, CMD
mean_diffs2 <- sp2 %>% 
  group_by(Trial,Provenance) %>% 
  summarise(meanH = mean(H, na.rm=TRUE),
            MAP = mean(MAPdiff),
            MSP = mean(MSPdiff),
            DD0 = mean(DD0diff),
            DD18 = mean(DD18diff),
            NFFD = mean(DD18diff),
            bFFP = mean(bFFPdiff),
            FFP = mean(FFPdiff),
            eFFP = mean(eFFPdiff),
            PAS = mean(PASdiff),
            Eref = mean(Erefdiff),
            CMD = mean(CMDdiff))

mean_diffs2 %>% 
  gather(key='variable',value='mean.diff', -Trial,-Provenance,-meanH) %>% 
  ggplot(aes(mean.diff,meanH))+
  facet_wrap(~variable)+
  geom_point()+
  stat_smooth()+
  ylab("Mean height (mm)")+
  xlab("Mean difference between trial and provenance sites")+
  ggtitle("Height ~ variables with large variation between trial and prov sites")

# re-do with standardised vars
sp3 <- read.csv("./Scots_pine/Scots_pine_H_cent_scal_allvars.csv") # 
sp3$X <- NULL
str(sp3)
colnames(sp3)[1]<-"H"
sp<-na.omit(sp)
sp3$Hraw <- sp$H
sp3$Trial <- sp$Trial
sp3$Provenance <- sp$Provenance
sp3 <- sp3 %>% mutate(MATdiff = MAT_P-MAT_T,
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
mean_diffs3 <- sp3 %>% 
  group_by(Trial,Provenance) %>% 
  summarise(meanH = mean(Hraw, na.rm=TRUE),
            MAP = mean(MAPdiff),
            MSP = mean(MSPdiff),
            DD0 = mean(DD0diff),
            DD18 = mean(DD18diff),
            NFFD = mean(DD18diff),
            bFFP = mean(bFFPdiff),
            FFP = mean(FFPdiff),
            eFFP = mean(eFFPdiff),
            PAS = mean(PASdiff),
            Eref = mean(Erefdiff),
            CMD = mean(CMDdiff))

mean_diffs3 %>% 
  gather(key='variable',value='mean.diff', -Trial,-Provenance,-meanH) %>% 
  ggplot(aes(mean.diff,meanH))+
  facet_wrap(~variable)+
  geom_point()+
  stat_smooth()+
  ylab("Mean height (mm)")+
  xlab("Mean difference between trial and provenance sites")+
  ggtitle("Height ~ standardised variables with large variation between trial and prov sites")

### similar to plot richard produced? use his code to reproduct plot with my climate data
psySum = psy %>%
  tidyr::pivot_longer(cols = x:VAPPressure)%>%
  group_by(Population, PlantingSite, name, value)%>%
  summarise(mn_ht = mean(na.omit(W17Height)))

ggplot(psySum, aes(value, mn_ht, colour = PlantingSite))+
  geom_point()+geom_smooth()+facet_wrap(~name, scales = "free")+
  theme_bw()+theme(legend.position = "bottom")+
  geom_smooth(data = psySum, aes(x = value, y = mn_ht, group = 1), lty = "dashed", method = "lm", colour = "black", se = F)

### lattice 
library(lattice)
histogram(~H,data=sp)
md3<- mean_diffs2 %>% 
  gather(key='variable',value='mean.diff', -Trial,-Provenance,-meanH) 
xyplot(meanH~mean.diff | variable, data=md3)


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
# calculate mean square errors (MSE)
mean(basic_model_differences ^ 2)
mean(mixed_model_differences ^ 2) # better predictions
# but using MSE to decide whether to include an explanatory variable in a linear model architecture
# isn't yet complete because of a problem: 
# whenever you use it you will find that the model with the additional explanatory variable 
# has smaller prediction errors than the base model! 
# the technique always gives the same indication: include the additional explanatory variable
