
# date: 17/02/20
# author: VB
# description: using Fagus climate and trait data, learn how to check for colinearity and choose the best variables for the model
# will need to do this for poplar and scots pine
# see what i can do with Fagus and compare to Marta's out the box model in Height_2070_RCP8.5.R

# tutorials/advice
# http://www.sthda.com/english/articles/39-regression-model-diagnostics/160-multicollinearity-essentials-and-vif-in-r/

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

# look at data

# correlation matrix
str(Fagus)
Fagus.cor <- cor(Fagus[,-c(1,3)], method = c("spearman")) # all variables expect factors
Fagus.cor<-as.data.frame(Fagus.cor)
#corrplot(Fagus.cor)

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
# Error in vif.default(model1) : there are aliased coefficients in the model
# Use the 'alias' function in R to see which variables are linearly dependent. Remove the dependent variables and the vif function should work correctly. https://stackoverflow.com/questions/28885160/vifs-returning-aliased-coefficients-in-r

# the linearly dependent variables
ld.vars <- attributes(alias(model1)$Complete)$dimnames[[1]]

# remove the linearly dependent variables 
Fagus<-Fagus[ , !(names(Fagus) %in% ld.vars)]
dim(Fagus)

training.samples <- Fagus$H %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- Fagus[training.samples, ]
test.data <- Fagus[-training.samples, ]

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

# remove random effects for now
fagus<-Fagus[,-c(1:7)] 
str(fagus)
training.samples <- fagus$H %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- fagus[training.samples, ]
test.data <- fagus[-training.samples, ]

model3 <- lm(H ~., data = train.data)
# Make predictions
predictions <- model3 %>% predict(test.data)
# Model performance
data.frame(
  RMSE = RMSE(predictions, test.data$H),
  R2 = R2(predictions, test.data$H)
)

# detect multicollinearity
VIF <- car::vif(model3) %>%
  as.list() %>% 
  as.data.frame() %>% 
  gather(key = 'variable', value = 'VIF') %>% 
  arrange(desc(VIF))

# remove variables with high VIF (above 5-10)
best.vars<-unique(VIF$variable[which(VIF$VIF<50)])
best.vars

# filter to just vars with low(ish) VIF
fagus.mod<-Fagus[,c(1:8,18,11,19,16,21,12)]
str(fagus.mod)

#marta_mod <- lmer(log(H) ~ St_Age + St_Pet.max_P + St_Bio13_T + I(St_Pet.max_P^2) + I(St_Bio13_T^2) +
             #St_Age*St_Pet.max_P + St_Age*St_Bio13_T + St_Pet.max_P*St_Bio13_T + 
             #(1|Trial/Block/Tree_ID) + (1|ID_ProvCode), Fagus)

fm1 <- lmer(H ~ pet.min_T + (1|Trial/Block/Tree_ID) + (1|ID_ProvCode), fagus.mod) # trial site min PET
summary(fm1)
par(mfrow = c(2,2))
plot(fm1)


# use coding club tutorial, replacing dragons dataset with fagus
# https://ourcodingclub.github.io/2017/03/15/mixed-models.html
# look at distribution of response variable
hist(Fagus$H)
# it is good practice to standardise your explanatory variables before proceeding...
# so that they have a mean of zero (“centering”) and standard deviation of one (“scaling”) 
# It ensures that the estimated coefficients are all on the same scale, making it easier to compare effect sizes.

Fagus$bio14_T_2 <- scale(Fagus$bio14_T, center = TRUE, scale = TRUE) # variable with lowest VIF, precipitation of driest month?
hist(Fagus$bio14_T_2)
basic.lm <- lm(H ~ bio14_T_2, data = Fagus)
summary(basic.lm)
# plot the data
(prelim_plot <- ggplot(Fagus, aes(x = bio14_T_2, y = H)) +
    geom_jitter() +
    geom_smooth(method = "lm"))
# plot residuals
plot(basic.lm, which = 1) # red line should be flat like the dashed grey line
# look at qqplot
plot(basic.lm, which = 2) # points should fall on diagonal line

# check for data independence - essentially checking if need to account for random effects
boxplot(H ~ Trial, data = Fagus)
(colour_plot <- ggplot(Fagus, aes(x = bio14_T_2, y = H, colour = Trial)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position = "none"))

(split_plot <- ggplot(aes(bio14_T_2, H), data = Fagus) + 
    geom_point() + 
    facet_wrap(~ Trial) + # create a facet for each mountain range
    xlab("precipitation driest month") + 
    ylab("height"))

# include trial as a fixed effect
trial.lm <- lm(H ~ bio14_T_2 + Trial, data = Fagus)
summary(trial.lm)

# as a random effect
mixed.lmer <- lmer(H ~ bio14_T_2 + (1|Trial), data = Fagus)
summary(mixed.lmer)
plot(mixed.lmer) 
qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer)) # points should fall on line

# account for nesting
mixed.lmer2 <- lmer(H ~ bio14_T_2 + (1|Trial/Block/Tree_ID), data=Fagus)
summary(mixed.lmer2)
plot(mixed.lmer2) 
qqnorm(resid(mixed.lmer2))
qqline(resid(mixed.lmer2)) # points should fall on line

# 22337/(22337 + 604 + 9292) # ~69 %
# trials explain 69% of the variation
# so the differences between trials explain ~69% of the variance that’s “left over” after the variance explained by bio14.

# account for nesting and provenance
mixed.lmer3 <- lmer(H ~ bio14_T_2 + (1|Trial/Block/Tree_ID) + (1|ID_ProvCode), data=Fagus)
summary(mixed.lmer3)
plot(mixed.lmer3) 
qqnorm(resid(mixed.lmer3))
qqline(resid(mixed.lmer3)) # points should fall on line - looks better. 

(mm_plot <- ggplot(Fagus, aes(x = bio14_T_2, y = H, colour = Trial)) +
    facet_wrap(~Trial, nrow=2) +   # a panel for each trial
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(Fagus, pred = predict(mixed.lmer3)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)


# some useful code from datacamp
tidy(mixed.lmer3)
# save the model predictions as a column to the original data.frame
Fagus.s$lmerPredict <- predict(mixed.lmer3)
# plot the original data
ggFagus2 <- ggplot( Fagus.s, aes(x = St_Bio14_P, y = H) ) +
  geom_point()+
  theme_minimal() +
  geom_abline(data = Fagus.s,
              aes(intercept = intercept, slope = slope))
# use the predicted values to plot the new outputs
ggFagus2 +
  geom_line( data =  Fagus,
             aes(x = x, y = lmerPredict, color = Trial),
             linetype = 2)
