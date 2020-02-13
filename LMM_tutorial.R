
# https://www.r-bloggers.com/linear-mixed-effect-models-in-r/

# Install (if necessary) and load nlme and lme4
library(nlme)
library(lme4)

# load data
Fagus <- read.csv("~/R/DeltaTraitSDM/Fagus/Fagus_DeltaTraitSDM_H.csv",header=T)
head(Fagus)
str(Fagus)
dim(Fagus) # 203105 observations, 77 variables

# correlation matrix
Fagus.s <- Fagus[,-c(15:77)]
head(Fagus.s)
Fagus.s <- na.omit(Fagus.s)
str(Fagus.s)
Fagus.s <- Fagus.s[,-c(1:11)]
Fagus.cor <- cor(Fagus.s[,-c(1:3)], method = c("spearman"))
library(corrplot)
corrplot(Fagus.cor)


##################################################################################################
# detect multicollinearity using VIF http://www.sthda.com/english/articles/39-regression-model-diagnostics/160-multicollinearity-essentials-and-vif-in-r/
##################################################################################################

library(tidyverse)
library(caret)

# Split the data into training and test set
set.seed(123)
training.samples <- Fagus.s$H %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- Fagus.s[training.samples, ]
test.data <- Fagus.s[-training.samples, ]

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
# Error in vif.default(model1) : there are aliased coefficients in the model
# Use the 'alias' function in R to see which variables are linearly dependent. Remove the dependent variables and the vif function should work correctly. https://stackoverflow.com/questions/28885160/vifs-returning-aliased-coefficients-in-r

# the linearly dependent variables
ld.vars <- attributes(alias(model1)$Complete)$dimnames[[1]]

# remove the linearly dependent variables 
Fagus.s[ , !(names(Fagus.s) %in% ld.vars)]
dim(Fagus.s)

training.samples <- Fagus.s$H %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- Fagus.s[training.samples, ]
test.data <- Fagus.s[-training.samples, ]

model2 <- lm(H ~., data = train.data)
# Make predictions
predictions <- model2 %>% predict(test.data)
# Model performance
data.frame(
  RMSE = RMSE(predictions, test.data$H),
  R2 = R2(predictions, test.data$H)
)

# detect multicollinearity
car::vif(model2)
# same issue.

# make my own smaller dataset

fagusP<-Fagus[,c(1,12,13,17:37)] # just provenance data
training.samples <- fagusP$H %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- fagusP[training.samples, ]
test.data <- fagusP[-training.samples, ]

modelProv <- lm(H ~., data = train.data)
# Make predictions
predictions <- modelProv %>% predict(test.data)
# Model performance
data.frame(
  RMSE = RMSE(predictions, test.data$H),
  R2 = R2(predictions, test.data$H)
)

# detect multicollinearity
VIF_P<-car::vif(modelProv)
VIF_P<-as.data.frame(VIF_P)
arrange(desc(VIF_P$VIF_P))

# remove variables with high VIF (above 5-10)

##################################################################################################
# visualise distributions of variables 
##################################################################################################

# Overview of the variables
par(mfrow = c(2,4))
barplot(table(reg), ylab = "Frequency", main = "Region")
barplot(table(popu), ylab = "Frequency", main = "Population")
barplot(table(gen), ylab = "Frequency", las = 2, main = "Genotype")
barplot(table(rack), ylab = "Frequency", main = "Rack")
barplot(table(nutrient), ylab = "Frequency", main = "Nutrient")
barplot(table(amd), ylab = "Frequency", main = "AMD")
barplot(table(status), ylab = "Frequency", main = "Status")
hist(total.fruits, col = "grey", main = "Total fruits", xlab = NULL)