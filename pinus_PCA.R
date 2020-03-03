
# date: 03/03/20
# author: VB
# description: using Scots pine climate and trait data, carry out PCA to find best variables

# libraries
library(tidyverse)
library(car)
library(corrplot)
library(caret)
library(nlme)
library(lme4)
library(broom)
library(ggplot2)

# read in data
Pinus<-read.csv("./Scots_pine/Scots_pine_H.csv")
Pinus$X<-NULL
# remove NAs
Pinus<-na.omit(Pinus)
str(Pinus)

#####################################################
# PCA
#####################################################
Pinus<-Pinus[,-c(1:10,12)]
str(Pinus)
Pinus$DD18_T <- NULL
corrplot(cor(Pinus), method = "ellipse")
Pinus<-Pinus[complete.cases(Pinus),]
Pinus.pca <- prcomp(Pinus, cor=TRUE, scores=TRUE) # https://stackoverflow.com/questions/26396356/princomp-error-in-r-covariance-matrix-is-not-non-negative-definite
summary(Pinus.pca)

# choose components to use based on proportion of variance explained
screeplot(Pinus.pca, type = "lines")
# or by cumulative variance
# Variance explained
pc.var <- Pinus.pca$sdev^2
# Proportion of variation
pc.pvar <- pc.var / sum(pc.var)
# Cumulative proportion
plot(cumsum(pc.pvar), type = 'b')
abline(h = 0.9)

attributes(Pinus.pca)
loadings(Pinus.pca)

library(factoextra)
fviz_pca_ind(Pinus.pca)
fviz_pca_var(Pinus.pca)
fviz_pca_biplot(Pinus.pca)

# try post transformation
# transform variables with large values as they will dominate any correlation coefficient
Pinus$W17Height<-log(Pinus$W17Height)
Pinus$Elevation<-log(Pinus$Elevation)
Pinus$MAP_P<-log(Pinus$MAP_P)
Pinus$MSP_P<-log(Pinus$MSP_P)
Pinus$DD0_P<-log(Pinus$DD0_P)
Pinus$DD5_P<-log(Pinus$DD5_P)
Pinus$DD_18_P<-log(Pinus$DD_18_P)
Pinus$NFFD_P<-log(Pinus$NFFD_P)
Pinus$bFFP_P<-log(Pinus$bFFP_P)
Pinus$eFFP_P<-log(Pinus$eFFP_P)
Pinus$FFP_P<-log(Pinus$FFP_P)
Pinus$PAS_P<-log(Pinus$PAS_P)
Pinus$Eref_P<-log(Pinus$Eref_P)
Pinus$MAP_T<-log(Pinus$MAP_T)
Pinus$MSP_T<-log(Pinus$MSP_T)
Pinus$DD0_T<-log(Pinus$DD0_T)
Pinus$DD5_T<-log(Pinus$DD5_T)
Pinus$DD_18_T<-log(Pinus$DD_18_T)
Pinus$NFFD_T<-log(Pinus$NFFD_T)
Pinus$bFFP_T<-log(Pinus$bFFP_T)
Pinus$eFFP_T<-log(Pinus$eFFP_T)
Pinus$FFP_T<-log(Pinus$FFP_T)
Pinus$PAS_T<-log(Pinus$PAS_T)
Pinus$Eref_T<-log(Pinus$Eref_T)
Pinus$CMD_P<-log(Pinus$CMD_P)
Pinus$CMD_P[which(Pinus$CMD_P==-Inf)]<-0
Pinus$CMD_T<-log(Pinus$CMD_T)
Pinus$CMD_T[which(Pinus$CMD_T==-Inf)]<-0

head(Pinus)
summary(Pinus)

Pinus.pca2 <- prcomp(Pinus, cor=TRUE, scores=TRUE) 
# choose components to use based on proportion of variance explained
screeplot(Pinus.pca2, type = "lines")
# or by cumulative variance
# Variance explained
pc.var <- Pinus.pca2$sdev^2
# Proportion of variation
pc.pvar <- pc.var / sum(pc.var)
# Cumulative proportion
plot(cumsum(pc.pvar), type = 'b')
abline(h = 0.9)

attributes(Pinus.pca2)
loadings(Pinus.pca2)
fviz_pca_var(Pinus.pca2)
