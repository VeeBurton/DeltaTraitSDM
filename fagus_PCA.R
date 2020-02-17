

# date: 17/02/20
# author: VB
# description: using Fagus climate and trait data, carry out PCA to find best variables

# libraries
library(tidyverse)
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
Fagus<-Fagus[,-c(11:13)] # lots of NAs
Fagus<-na.omit(Fagus)

# detect multicollinearity using VIF 
# split the data into training and test set
training.samples <- Fagus$H %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- Fagus[training.samples, ]
test.data <- Fagus[-training.samples, ]

# build a regression model with all variables
model1 <- lm(H ~., data = train.data)
# make predictions
predictions <- model1 %>% predict(test.data)
# model performance
data.frame(
  RMSE = RMSE(predictions, test.data$H),
  R2 = R2(predictions, test.data$H)
)

# detect multicollinearity
car::vif(model1)
# the linearly dependent variables
ld.vars <- attributes(alias(model1)$Complete)$dimnames[[1]]

#####################################################
# variance and covariance
#####################################################

str(Fagus)
Fagus$Year_measurement<-as.numeric(Fagus$Year_measurement)
Fagus$Tree_ID<-as.numeric(Fagus$Tree_ID)
Fagus$ProvCode<-as.numeric(Fagus$ProvCode)
Fagus$Block<-as.numeric(Fagus$Block)
Fagus$Tree<-as.numeric(Fagus$Tree)
Fagus$YearPlantation<-as.numeric(Fagus$YearPlantation)
Fagus$Age<-as.numeric(Fagus$Age)
#Fagus.s<-Fagus[,-c(1,3)]
str(Fagus)

library(corrplot)
var(Fagus.s)
cor(Fagus.s)
corrplot(cor(Fagus.s), method = "ellipse")


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