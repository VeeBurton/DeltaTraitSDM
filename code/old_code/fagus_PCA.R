

# date: 17/02/20
# author: VB
# description: using Fagus climate and trait data, carry out PCA to find best variables

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
Fagus <- read.table("./Fagus/Fagus_H_Final.txt",header=T)
str(Fagus)
Fagus<-Fagus[,-c(16:72)]
Fagus<-Fagus[,-c(8:9)]
summary(Fagus)
Fagus<-Fagus[,-c(11:13)] # lots of NAs - traits DBH, BasalD, Mort
Fagus<-na.omit(Fagus)


#####################################################
# PCA
#####################################################
Fagus.s<-Fagus.s[,-c(11:55)]
corrplot(cor(Fagus.s), method = "ellipse")
Fagus.s<-Fagus.s[complete.cases(Fagus.s),]
fagus.pca <- prcomp(Fagus.s, cor=TRUE, scores=TRUE) # https://stackoverflow.com/questions/26396356/princomp-error-in-r-covariance-matrix-is-not-non-negative-definite
summary(fagus.pca)

# choose components to use based on proportion of variance explained
screeplot(fagus.pca, type = "lines")
# or by cumulative variance
# Variance explained
pc.var <- fagus.pca$sdev^2
# Proportion of variation
pc.pvar <- pc.var / sum(pc.var)
# Cumulative proportion
plot(cumsum(pc.pvar), type = 'b')
abline(h = 0.9)

attributes(fagus.pca)
loadings(fagus.pca)

library(factoextra)
fviz_pca_ind(fagus.pca)
fviz_pca_var(fagus.pca)
fviz_pca_biplot(fagus.pca)