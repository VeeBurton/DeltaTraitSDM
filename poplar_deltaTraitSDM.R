
# libraries
library(tidyverse)
library(raster)
library(rgdal) 
library(rworldmap)
library(sjPlot) 
library(lattice)
library(spatial.tools) 
library(maptools)
library(dismo)

setwd("C:/Users/vanessa.burton/OneDrive - Forest Research/Documents/R/DeltaTraitSDM")

poplar<-read.csv('./Poplar/poplar_monthly_processed.csv')

head(poplar)
str(poplar)
dim(poplar) # 203105 observations, 77 variables

# run through biovars?
poplar.b<- poplar %>%
  group_by(site) %>%
  do(data.frame(biovars(.$pr, .$tasmin, .$tasmax)))

poplar.all<-merge(poplar,poplar.b,by='site')
str(poplar.all)

##############################################
# ALSO NEED TRAIT DATA TO MAKE THIS WORTHWHILE
##############################################

# correlation matrix
library(corrplot)
poplar.s <- poplar.all[,-c(1,2,3,7)]
head(poplar.s)
poplar.s <- na.omit(poplar.s)
str(poplar.s)
poplar.cor <- cor(poplar.s, method = c("spearman"))
corrplot(poplar.cor)

# detect multicollinearity using VIF
library(caret)

# Split the data into training and test set
set.seed(123)
training.samples <- poplar.s$H %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- poplar.s[training.samples, ]
test.data <- poplar.s[-training.samples, ]

# Build a regression model with all variables
model1 <- lm(H ~., data = train.data)
# Make predictions
predictions <- model1 %>% predict(test.data)
# Model performance
data.frame(
  RMSE = RMSE(predictions, test.data$H),
  R2 = R2(predictions, test.data$H)
)

# detect multicollinearity
car::vif(model1)
