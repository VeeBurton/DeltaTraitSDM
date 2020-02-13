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
library(tidyverse)
library(broom)
###################

#################################################################################################
# PLAN: scots pine
# test height against a number of most likely explanatory variables:
# growing season length?
# frost-free period?
# start with these and see if relationships are strong enough to build a more complex model?

sp <- read.csv()

# plot data by trial/provenance etc.
ggplot(aes(x, height), data = sp) + geom_point() +
  facet_wrap(~ x) +
  xlab("x") + ylab("height")

boxplot(height ~ x, data = sp)

# Choosing variables
# Stepwise modelling
# Combination of variables (PCA etc)
# Look at variable colinearity

# look at distribution of height (response variable)
hist(sp$height) 

# good practice to standardise all explanatory variables (x)
# include Trial/Block/Tree_id and Provenance as random effects - (1|Trial/Block/Tree_id)
# may need to make nesting explicit by creating new variables e.g. Trial1-Blocka, Trial1-Blockb etc.

# e.g.
SPmod1 <- lmer(height ~ growing.season 
               + FFP 
               + (1|Trial/Block/Tree_id))

summary(SPmod1)
coef(Spmod1)
augment(Spmod1)

# plot residuals
plot(SPmod1, which = 1)
# q plot
qqnorm(resid(SPmod1))
qqline(resid(SPmod1))

