
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
library(factoextra)

# read in data
Pinus<-read.csv("./Scots_pine/Scots_pine_H.csv")
Pinus$X<-NULL
# remove NAs
Pinus<-na.omit(Pinus)
str(Pinus)

#####################################################
# PCA (on original data)
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

#fviz_pca_ind(Pinus.pca)
fviz_pca_var(Pinus.pca, repel=TRUE)
#fviz_pca_biplot(Pinus.pca, repel = TRUE)

#####################################################
# try post transformation
# log only variables with large values as they will dominate any correlation coefficient
#####################################################
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

#####################################################
# log all vars
#####################################################
Pinus<-read.csv("./Scots_pine/Scots_pine_H.csv")
Pinus$X<-NULL
# remove NAs
Pinus<-na.omit(Pinus)
str(Pinus)
Pinus<-Pinus[,-c(1:10,12)]
variables<-unique(colnames(Pinus))

for (i in c(1:44)){
  #i<-1
  var<-variables[i]
  var1<-Pinus[, c(var)] 
  var1<-as.numeric(var1)
  Pinus[,c(var)]<-log(var1)
}

summary(Pinus)

Pinus$Longitude<-NULL
Pinus$DD18_P<-NULL
Pinus$EMNT_P<-NULL
Pinus$DD18_T<-NULL
Pinus$EMNT_T<-NULL
summary(Pinus)
Pinus$CMD_T[which(Pinus$CMD_T==-Inf)]<-0
Pinus$CMD_P[which(Pinus$CMD_P==-Inf)]<-0
summary(Pinus)
Pinus<-na.omit(Pinus)

Pinus.pca3 <- prcomp(Pinus, cor=TRUE, scores=TRUE) 
summary(Pinus.pca3)
# choose components to use based on proportion of variance explained
screeplot(Pinus.pca3, type = "lines")
# or by cumulative variance
# Variance explained
pc.var <- Pinus.pca3$sdev^2
# Proportion of variation
pc.pvar <- pc.var / sum(pc.var)
# Cumulative proportion
plot(cumsum(pc.pvar), type = 'b')
abline(h = 0.9)

attributes(Pinus.pca3)
loadings(Pinus.pca3)
fviz_pca_var(Pinus.pca3)

#####################################################
# standardise all
#####################################################

Pinus<-read.csv("./Scots_pine/Scots_pine_H.csv")
Pinus$X<-NULL
# remove NAs
Pinus<-na.omit(Pinus)
str(Pinus)
Pinus<-Pinus[,-c(1:10,12)]
variables<-unique(colnames(Pinus))

for (i in c(1:44)){
  #i<-1
  var<-variables[i]
  var1<-Pinus[, c(var)] 
  var1<-as.numeric(var1)
  Pinus[,c(var)]<-robustHD::standardize(var1, centerFun = mean)
}

summary(Pinus)
Pinus$DD18_T<-NULL

Pinus.pca4 <- prcomp(Pinus, cor=TRUE, scores=TRUE) 
summary(Pinus.pca4)
# choose components to use based on proportion of variance explained
screeplot(Pinus.pca4, type = "lines")
# or by cumulative variance
# Variance explained
pc.var <- Pinus.pca4$sdev^2
# Proportion of variation
pc.pvar <- pc.var / sum(pc.var)
# Cumulative proportion
plot(cumsum(pc.pvar), type = 'b')
abline(h = 0.9)

attributes(Pinus.pca4)
loadings(Pinus.pca4)
fviz_pca_var(Pinus.pca4, repel = TRUE)

write.csv(Pinus, "./Scots_pine/Scots_pine_H_standardised_allvars.csv" )

###############################
# centred and scaled
###############################
Pinus<-read.csv("./Scots_pine/Scots_pine_H.csv")
Pinus$X<-NULL

# remove NAs
Pinus<-na.omit(Pinus)
head(Pinus)
Pinus<-Pinus[,-c(1:10,12)]
variables<-unique(colnames(Pinus))

c. <- function (x) scale(x, scale = FALSE)

z. <- function (x) scale(x)

# using functions from diagnosing_collinearity.R 
for (i in c(1:44)){
  #i<-1
  var1<-variables[i]
  var2<-Pinus[, c(var1)] 
  var_cent<-c.(var2)
  var_scale<-z.(var_cent)
  Pinus[,c(var1)]<-var_scale
}

summary(Pinus)
Pinus$DD18_T.V1<-NULL

Pinus.pca5 <- prcomp(Pinus, cor=TRUE, scores=TRUE) 
summary(Pinus.pca5)
# choose components to use based on proportion of variance explained
screeplot(Pinus.pca5, type = "lines")
# or by cumulative variance
# Variance explained
pc.var <- Pinus.pca5$sdev^2
# Proportion of variation
pc.pvar <- pc.var / sum(pc.var)
# Cumulative proportion
plot(cumsum(pc.pvar), type = 'b')
abline(h = 0.9)

attributes(Pinus.pca5)
loadings(Pinus.pca5)
fviz_pca_var(Pinus.pca5, repel = TRUE)

write.csv(Pinus, "./Scots_pine/Scots_pine_H_cent_scal_allvars.csv" )

# remove lat/long/elev to see if that changes things
head(Pinus)
Pinus<-Pinus[,-c(2:4)]
Pinus.pca6 <- prcomp(Pinus, cor=TRUE, scores=TRUE)
screeplot(Pinus.pca6, type = "lines")
fviz_pca_var(Pinus.pca6, repel = TRUE)
# no major changes

# more PCA analysis (of centred and scaled variables)


