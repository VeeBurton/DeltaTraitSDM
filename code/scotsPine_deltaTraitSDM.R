# libraries
###################
library(lme4)
library(nlme)
library(MuMIn)
library(raster)
library(rgdal)
library(rworldmap)
library(sjPlot) 
library(curry)
library(broom)
library(tidyverse)
library(corrplot)
library(caret)
library(ggplot2)
library(performance)
###################

# based on all analysis so far, these are the most important variables according to:

# pairplots:
    # MWMT_T (0.41)
    # PAS_T (0.41)
    # TD_T (0.33)
    # DD0_T (0.30)
    # MCMT_T (0.28)
    # eFFP_T (0.28)
    # Eref_T (0.24)
    # EMT_T (0.22)
    # DD18_T (0.21)
    # NFFD_T (0.21) 
# PCA: 
    # DD_18_T
    # FFP_T
    # bFFP_P
    # FFP_P
    # TD_P
# exploration/difference between trial & provenance:
    # PAS
    # TD
    # MWMT
    # DD0
    # MCMT

# read in standardised data (with raw height)
sp.raw <- read.csv("./Scots_pine/Scots_pine_H.csv") # raw data
sp.raw$X <- NULL
str(sp.raw)
sp.raw<-na.omit(sp.raw)
sp <- read.csv("./Scots_pine/Scots_pine_H_cent_scal_allvars.csv") # 
sp$X <- NULL
str(sp)
colnames(sp)[1]<-"Hscal"
sp$Hraw <- sp.raw$W17Height
sp$Trial <- sp.raw$Trial
sp$Provenance <- sp.raw$Provenance
sp$Block <- sp.raw$Block
colnames(sp)[22]<-'EMT_P'
colnames(sp)[41]<-'EMT_T'
rm(sp.raw)

# models to compare
pair.mod1 <- lmer(Hraw ~ MWMT_T + PAS_T + (1|Trial/Block/Provenance), data=sp)
pair.mod2 <- lmer(log(Hraw) ~ MWMT_T + PAS_T + (1|Trial/Block/Provenance), data=sp)
pair.mod3 <- lmer(log(Hraw) ~ TD_T + Eref_T + (1|Trial/Block/Provenance), data=sp)
pca.mod1 <- lmer(log(Hraw) ~ DD_18_T + FFP_P + (1|Trial/Block/Provenance), data=sp)
pca.mod2 <- lmer(log(Hraw) ~ DD_18_T + TD_P + (1|Trial/Block/Provenance), data=sp)
pca.mod3 <- lmer(log(Hraw) ~ FFP_T + TD_P + (1|Trial/Block/Provenance), data=sp)
diff.mod1 <- lmer(log(Hraw) ~ PAS_P + TD_T + (1|Trial/Block/Provenance), data=sp)
diff.mod2 <- lmer(log(Hraw) ~ MWMT_P + MCMT_T + (1|Trial/Block/Provenance), data=sp)
diff.mod3 <- lmer(log(Hraw) ~ DD0_P + PAS_T + (1|Trial/Block/Provenance), data=sp)

# performance package
# comparisons
compare_performance(pair.mod1,pair.mod2,pair.mod3,pca.mod1,pca.mod2,pca.mod3,diff.mod1,diff.mod2,diff.mod3, rank = TRUE)
#Model     |    Type |      AIC |      BIC | R2_marginal |  RMSE | Performance_Score
#-----------------------------------------------------------------------------------
#diff.mod3 | lmerMod |  4652.06 |  4690.49 |        0.17 |  0.88 |            99.99%
#pair.mod3 | lmerMod |  4651.48 |  4689.92 |        0.17 |  0.87 |            99.82%
#pair.mod2 | lmerMod |  4650.35 |  4688.78 |        0.17 |  0.88 |            99.69%
#diff.mod1 | lmerMod |  4659.38 |  4697.82 |        0.10 |  0.87 |            87.85%
#diff.mod2 | lmerMod |  4658.95 |  4697.39 |        0.07 |  0.87 |            83.19%
#pca.mod2  | lmerMod |  4657.62 |  4696.06 |        0.05 |  0.87 |            78.22%
#pca.mod1  | lmerMod |  4654.09 |  4692.52 |        0.03 |  0.88 |            76.19%
#pca.mod3  | lmerMod |  4657.69 |  4696.12 |        0.03 |  0.87 |            74.98%
#pair.mod1 | lmerMod | 21221.14 | 21259.58 |        0.12 | 90.20 |            15.98%
plot(compare_performance(pca.mod1,pca.mod2,pca.mod3,diff.mod3, rank = TRUE))
r2(diff.mod3)$R2_marginal
icc(diff.mod3)
check_model(diff.mod3)

# "best models" from dredge
dredge1<-lmer(log(Hraw) ~ DD_18_T + NFFD_T + (1|Trial/Block/Provenance),data = sp)
dredge2<-lmer(log(Hraw) ~ DD_18_T + NFFD_T + PAS_T + (1|Trial/Block/Provenance),data = sp)
dredge3<-lmer(log(Hraw) ~ DD_18_T + NFFD_T + TD_T + (1|Trial/Block/Provenance),data = sp)

compare_performance(diff.mod3,pair.mod3,pair.mod2,dredge1,dredge2,dredge3, rank = TRUE)
#Model     |    Type |     AIC |     BIC | R2_marginal | RMSE |    BF | Performance_Score
#----------------------------------------------------------------------------------------
#dredge1   | lmerMod | 4646.35 | 4684.78 |        0.17 | 0.88 | 17.37 |            60.00%
#dredge2   | lmerMod | 4646.35 | 4684.78 |        0.17 | 0.88 | 17.37 |            60.00%
#dredge3   | lmerMod | 4646.35 | 4684.78 |        0.17 | 0.88 | 17.37 |            60.00%
#diff.mod3 | lmerMod | 4652.06 | 4690.49 |        0.17 | 0.88 |  1.00 |            37.77%
#pair.mod3 | lmerMod | 4651.48 | 4689.92 |        0.17 | 0.87 |  1.33 |            33.33%
#pair.mod2 | lmerMod | 4650.35 | 4688.78 |        0.17 | 0.88 |  2.35 |            13.64%
plot(compare_performance(diff.mod3,pair.mod3,pair.mod2,dredge1,dredge2,dredge3, rank = TRUE))

whittet.mod1 <- lmer(log(Hraw) ~ MAP_T + FFP_P + (1|Trial/Block/Provenance), data = sp)
whittet.mod2 <- lmer(log(Hraw) ~ MAP_P + FFP_T  + (1|Trial/Block/Provenance), data = sp)
whittet.mod3 <- lmer(log(Hraw) ~ DD5_T + MWMT_P + MCMT_T + MAP_T + TD_P + (1|Trial/Block/Provenance), data = sp)
compare_performance(dredge1,diff.mod3,pair.mod3,pair.mod2,whittet.mod1,whittet.mod2,whittet.mod3, rank=TRUE)
#Model        |    Type |     AIC |     BIC | R2_marginal | RMSE |   BF | Performance_Score
#------------------------------------------------------------------------------------------
#dredge1      | lmerMod | 4646.35 | 4684.78 |        0.17 | 0.88 | 1.00 |            79.79%
#pair.mod3    | lmerMod | 4651.48 | 4689.92 |        0.17 | 0.87 | 0.08 |            71.42%
#diff.mod3    | lmerMod | 4652.06 | 4690.49 |        0.17 | 0.88 | 0.06 |            67.81%
#pair.mod2    | lmerMod | 4650.35 | 4688.78 |        0.17 | 0.88 | 0.14 |            54.72%
#whittet.mod2 | lmerMod | 4653.67 | 4692.11 |        0.04 | 0.88 | 0.03 |            45.52%
#whittet.mod1 | lmerMod | 4654.40 | 4692.84 |        0.00 | 0.88 | 0.02 |            40.66%
#whittet.mod3 | lmerMod | 4662.79 | 4712.21 |        0.10 | 0.87 | 0.00 |            30.82%
plot(compare_performance(dredge1,diff.mod3,pair.mod3,pair.mod2,whittet.mod1,whittet.mod2,whittet.mod3,rank=TRUE))

# run on training and test datasets
results<-NULL

for(i in 1:1000) {
  
  #i<-1
  rows<-sample(x=1:nrow(sp),size=round(nrow(sp)*0.25,1))
  training<-sp[-rows,]
  testing<-sp[rows,]
  
  dredge1<-lmer(log(Hraw) ~ DD_18_T + NFFD_T + (1|Trial/Block/Provenance),data = training)
  pair.mod2 <- lmer(log(Hraw) ~ MWMT_T + PAS_T + (1|Trial/Block/Provenance), data=training)
  diff.mod3 <- lmer(log(Hraw) ~ DD0_P + PAS_T + (1|Trial/Block/Provenance), data=training)
  pair.mod2 <- lmer(log(Hraw) ~ MWMT_T + PAS_T + (1|Trial/Block/Provenance), data=training)
  whittet.mod2 <- lmer(log(Hraw) ~ MAP_P + FFP_T  + (1|Trial/Block/Provenance), data = training)
  
  p.dr1<-predict(dredge1,testing)-log(testing$Hraw)
  p.pr2<-predict(pair.mod2,testing)-log(testing$Hraw)
  p.df3<-predict(diff.mod3,testing)-log(testing$Hraw)
  p.pr2<-predict(pair.mod2,testing)-log(testing$Hraw)
  p.wh2<-predict(whittet.mod2,testing)-log(testing$Hraw)
  
  results<-rbind(results,data.frame(Run=rep(i,5),Model=c('dredge1','pair.mod2','diff.mod3','pair.mod2','whittet.mod2'),
                                    MAE=c(mean(abs(p.dr1)),mean(abs(p.pr2)),mean(abs(p.df3)),mean(abs(p.pr2)),mean(abs(p.wh2))),
                                    R2marg=c(r2(dredge1)$R2_marginal, r2(pair.mod2)$R2_marginal, r2(diff.mod3)$R2_marginal,
                                             r2(pair.mod2)$R2_marginal,r2(whittet.mod2)$R2_marginal)))
  
}

summary(results)
hist(results$MAE,breaks=50)
hist(results$R2marg,breaks=50)

results %>% 
  ggplot(aes(Model,MAE, colour=Model))+geom_boxplot()+
  theme_minimal()+
  ylab("MAE")+
  theme(axis.text.x = element_blank())

results %>% 
  ggplot(aes(Model,R2marg, colour=Model))+geom_boxplot()+
  theme_minimal()+
  ylab("R2marg")+
  theme(axis.text.x = element_blank())

resAVE<-aggregate(MAE~Model,results,FUN=mean)
resAVE<-resAVE[order(resAVE[,2],decreasing=F),]

