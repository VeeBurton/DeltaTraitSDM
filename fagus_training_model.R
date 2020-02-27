
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
head(Fagus)
# just random effects and standardised variables
Fagus<-Fagus[,c(1:7,12,119:135)]
summary(Fagus)
Fagus<-na.omit(Fagus)

# look at data

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
model1_summary<-tidy(model1)

# detect multicollinearity
car::vif(model1) # there are aliased coefficients in the model
summary(model1)$coeff
length(unique(summary(model1)$coeff))
alias(model1)$Complete
# the linearly dependent variables
ld.vars <- attributes(alias(model1)$Complete)$dimnames[[1]]
ld.vars
#Fagus$Year_measurement<-NULL
#Fagus$YearPlantation<-NULL
#Fagus$Trial<-NULL
#Fagus$Tree_ID<-NULL
#Fagus$ProvCode<-NULL

#####################################################
# variance and covariance
#####################################################

str(Fagus)
Fagus$Year_measurement<-as.numeric(Fagus$Year_measurement)
Fagus$Tree_ID<-as.numeric(Fagus$Tree_ID)
Fagus$ProvCode<-as.numeric(Fagus$ProvCode)
Fagus$Block<-as.numeric(Fagus$Block)
Fagus$Tree<-as.numeric(Fagus$Tree)
#Fagus$YearPlantation<-as.numeric(Fagus$YearPlantation)
Fagus$St_Age<-as.numeric(Fagus$St_Age)
str(Fagus)

Fagus.s<-Fagus[,-c(1,3)] # non-numeric
var(Fagus.s)
cor(Fagus.s)
corrplot(cor(Fagus.s), method = "ellipse")

training.samples <- Fagus.s$H %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- Fagus.s[training.samples, ]
test.data <- Fagus.s[-training.samples, ]

# build a regression model with all variables
model2 <- lm(H ~., data = train.data)
# make predictions
predictions <- model2 %>% predict(test.data)
# model performance
data.frame(
  RMSE = RMSE(predictions, test.data$H),
  R2 = R2(predictions, test.data$H)
)
model2_summary<-tidy(model2)

# detect multicollinearity
car::vif(model2) 
summary(model2)$coeff

# detect multicollinearity
VIF <- car::vif(model2) %>%
  as.list() %>% 
  as.data.frame() %>% 
  gather(key = 'variable', value = 'VIF') %>% 
  arrange(desc(VIF))

# remove variables with high VIF (above 5-10)
best.vars<-unique(VIF$variable[which(VIF$VIF<12)])
best.vars

# filter to just vars with low VIF
fagus.mod<-Fagus[,c('H', best.vars)]
str(fagus.mod)

#marta_mod <- lmer(log(H) ~ St_Age + St_Pet.max_P + St_Bio13_T + I(St_Pet.max_P^2) + I(St_Bio13_T^2) +
             #St_Age*St_Pet.max_P + St_Age*St_Bio13_T + St_Pet.max_P*St_Bio13_T + 
             #(1|Trial/Block/Tree_ID) + (1|ID_ProvCode), Fagus)

# use coding club tutorial, replacing dragons dataset with fagus
# https://ourcodingclub.github.io/2017/03/15/mixed-models.html
# look at distribution of response variable
hist(fagus.mod$H)
# it is good practice to standardise your explanatory variables before proceeding...
# so that they have a mean of zero (“centering”) and standard deviation of one (“scaling”) 
# It ensures that the estimated coefficients are all on the same scale, making it easier to compare effect sizes.
fagus.mod$St_H <- scale(fagus.mod$H, center = TRUE, scale = TRUE) # variable with lowest VIF, precipitation of driest month?
hist(fagus.mod$St_H)
# or log
log(fagus.mod$H)

basic.lm <- lm(H ~ St_Age, data = fagus.mod)
summary(basic.lm)
# plot the data
(prelim_plot <- ggplot(fagus.mod, aes(x =St_Age, y = H)) +
    geom_jitter() +
    geom_smooth(method = "lm"))
# plot residuals
plot(basic.lm, which = 1) # red line should be flat like the dashed grey line
# look at qqplot
plot(basic.lm, which = 2) # points should fall on diagonal line

log.lm <- lm(log(H)~St_Age,data=fagus.mod)
summary(log.lm)
plot(log.lm, which = 1)
plot(log.lm, which = 2)

# check for data independence - essentially checking if need to account for random effects
boxplot(H ~ ProvCode, data = fagus.mod)
(colour_plot <- ggplot(fagus.mod, aes(x = St_Age, y = H, colour = ProvCode)) +
    geom_point(size = 2) +
    theme_classic() +
    theme(legend.position = "none"))

(split_plot <- ggplot(aes(St_Age, H), data = fagus.mod) + 
    geom_point() + 
    facet_wrap(~ ProvCode) + # create a facet for each trial
    xlab("age") + 
    ylab("height"))

# include provenance as a fixed effect
prov.lm <- lm(H ~ St_Age + ProvCode, data = fagus.mod)
prov.lm_summary<-tidy(trial.lm)

# as a random effect
mixed.lmer <- lmer(H ~ St_Age + (1|Trial), data = fagus.mod)
mixed.lmer_summary<-tidy(mixed.lmer)
plot(mixed.lmer) 
qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer)) # points should fall on line

# account for nesting
mixed.lmer2 <- lmer(H ~ St_Age + (1|ProvCode/Block/Tree_ID), data=fagus.mod)
mixed.lmer2_summary<-tidy(mixed.lmer2)
plot(mixed.lmer2) 
qqnorm(resid(mixed.lmer2))
qqline(resid(mixed.lmer2)) # points should fall on line

summary(mixed.lmer2)
# 3866/(3866 + 3149 + 482 + 2473) # ~38 %
# trials explain 38% of the variation
# so the differences between trials explain ~38% of the variance that’s “left over” after the variance explained by Age.

# account for nesting and provenance
mixed.lmer3 <- lmer(H ~ St_Age + (1|Trial/Block/Tree_ID) + (1|ID_ProvCode), data=fagus.mod)
mixed.lmer3_summary<-tidy(mixed.lmer3)
plot(mixed.lmer3) 
qqnorm(resid(mixed.lmer3))
qqline(resid(mixed.lmer3)) # points should fall on line 

(mm_plot <- ggplot(fagus.mod, aes(x = St_Age, y = H, colour = Trial)) +
    facet_wrap(~Trial) +   # a panel for each trial
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(fagus.mod, pred = predict(mixed.lmer3)), aes(y = pred), size = 1) +  # adding predicted line from mixed model 
    theme(legend.position = "none",
          panel.spacing = unit(2, "lines"))  # adding space between panels
)


# some useful code from datacamp
# tidy(mixed.lmer3)
# save the model predictions as a column to the original data.frame
fagus.mod$lmerPredict <- predict(mixed.lmer3)
# plot the original data
#ggFagus2 <- ggplot(fagus.mod, aes(x = St_Age, y = H)) +
  #geom_point()+
  #theme_minimal() +
  #geom_abline(data = fagus.mod,
              #aes(intercept = intercept, slope = slope))
# use the predicted values to plot the new outputs
#ggFagus2 +
  #geom_line( data = fagus.mod,
             #aes(x = x, y = lmerPredict, color = Trial),
             #linetype = 2)

library(ggeffects)

# Extract the prediction data frame
pred.mm <- ggpredict(mixed.lmer3, terms = c("St_Age"))  # this gives overall predictions for the model

# Plot the predictions 

(ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = fagus.mod,                      # adding the raw data (scaled values)
               aes(x = St_Age, y = H, colour = Trial)) + 
    labs(x = "Age", y = "Height", 
         title = "Age affects height across trials") + 
    theme_minimal()
)

ggpredict(mixed.lmer3, terms = c("St_Age", "Trial"), type = "re") %>% 
  plot() +
  labs(x = "Age", y = "Height", title = "Effect of age on height in Beech") + 
  theme_minimal()

library(stargazer)
stargazer(mixed.lmer3, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

library(broom) # also useful for looking at results

tidy()
augment(mixed.lmer3)
glance(mixed.lmer3)

td <- tidy(mixed.lmer3, conf.int = TRUE)
ggplot(td, aes(estimate, term, color = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high)) +
  geom_vline()

tidied_cv <- tidy(mixed.lmer3)
glance_cv <- glance(mixed.lmer3)
ggplot(tidied_cv, aes(lambda, estimate)) + geom_line(color = "red") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .2) +
  scale_x_log10() +
  geom_vline(xintercept = glance_cv$lambda.min) +
  geom_vline(xintercept = glance_cv$lambda.1se, lty = 2)

# model selection
# Focus on your question, don’t just plug in and drop variables from a model haphazardly until you make something “significant”. 
# Always choose variables based on biology/ecology
# Define your goals and questions and focus on that
# Also, don’t just put all possible variables in (i.e. don’t overfit)
# Remember that as a rule of thumb, you need 10 times more data than parameters you are trying to estimate

# Tests from worst to best:
# Wald Z-tests
# Wald t-tests (but LMMs need to be balanced and nested)
# Likelihood ratio tests (via anova() or drop1())
# MCMC or parametric bootstrap confidence intervals

# e.g. anova
anova(mixed.lmer2, mixed.lmer3)

# Entire model selection
# A few notes on the process of model selection. 
#There are two ways here: (i) “top-down”, where you start with a complex model and gradually reduce it
# (ii) “step up”, where you start with a simple model and add new variables to it. 
# Unfortunately, you might arrive at different final models by using those strategies and so you need to be careful.
# The model selection process recommended by Zuur et al. (2009) is a top-down strategy and goes as follows:
# fit a full model (he even recommends “beyond optimal” i.e. more complex than you’d expect or want it to be)
# sort out the random effects structure (use REML likelihoods or REML AIC or BIC)
# sort out fixed effects structure (either use REML the F-statistic or the t-statistic or compare nested ML models - keep your random effects constant)
# once you arrive at the final model present it using REML estimation

# Be mindful of what you are doing, prepare the data well and things should be alright


