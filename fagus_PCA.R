

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


#####################################################
# variance and covariance
#####################################################

str(nordic_clim)
nordic_clim$Elev<-as.numeric(nordic_clim$Elev)
nordic_clim<-na.omit(nordic_clim)

library(corrplot)
var(nordic_clim[,c(2,3,4,5,8,9,10,11,12,13)])
cor(nordic_clim[,c(2,3,4,5,8,9,10,11,12,13)])
corrplot(cor(nordic_clim[,c(2,3,4,5,8,9,10,11,12,13)]), method = "ellipse")

pairs(nordic_clim[,c(2,3,4,5,8,9,10,11,12,13)], col=nordic_clim$Country)

library(lattice)
splom(~ nordic_clim[,c(2,3,4,5,8,9,10,11,12,13)], col=nordic_clim$Country, pch=16)

library(GGally)
ggpairs(data=nordic_clim, columns = c(2,3,4,5,8,9,10,11,12,13), mapping = aes(color=Country))

library(scatterplot3d)
scatterplot3d(nordic_clim[,c(2,3,4,5,8,9,10,11,12,13)], color = as.numeric(nordic_clim$Country))

#####################################################
# PCA
#####################################################
nordic.sub<-nordic_clim[,-c(1,2,7,8,9)]
nordic.sub<-na.omit(nordic.sub)
nordic.pca <- princomp(nordic.sub, cor=TRUE, scores=TRUE)
summary(nordic.pca)

# choose components to use based on proportion of variance explained
screeplot(nordic.pca, type = "lines")
# or by cumulative variance
# Variance explained
pc.var <- nordic.pca$sdev^2
# Proportion of variation
pc.pvar <- pc.var / sum(pc.var)
# Cumulative proportion
plot(cumsum(pc.pvar), type = 'b')
abline(h = 0.9)

attributes(nordic.pca)
loadings(nordic.pca)

# plotting scores
scores <- data.frame(nordic.pca$scores)
ggplot(data = scores, aes(x = Comp.1, y = Comp.1, label = rownames(scores))) +
  geom_text(size = 4, col = "steelblue")

elev <- factor(nordic.sub$Elev)
ggplot(data = scores, aes(x = Comp.1, y = Comp.2, label = rownames(scores),
                          color = elev)) + geom_text(size = 1)+ theme(legend.position="none")

tn <- factor(nordic.sub$TN_jan_30yr)
ggplot(data = scores, aes(x = Comp.1, y = Comp.2, label = rownames(scores),
                          color = tn)) + geom_text(size = 1)+ theme(legend.position="none")

library(factoextra)
fviz_pca_ind(nordic.pca)
fviz_pca_var(nordic.pca)
fviz_pca_biplot(nordic.pca)