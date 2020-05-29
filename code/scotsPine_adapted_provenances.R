
# Objective
# show provenances and families better adapted to future abiotic and biotic stress

# wd
setwd("~/Documents/FR/FR_R/DeltaTraitSDM")

# libs
library(tidyverse)
library(ggplot2)
library(ggpubr) 
library(rstatix)
library(agricolae)
library(lme4)
library(broom)
library(performance)

# raw data
sp <- read.csv("data-raw/Scots_pine_H.csv")
glimpse(sp)
sp <- na.omit(sp)
sp$X<-NULL
colnames(sp)[11]<-"height"

# centre and scale climate variables
sp[16:55]<-scale(sp[16:55],center = TRUE, scale = TRUE)
glimpse(sp)

# summary stats
sp.prov <- sp %>%
  group_by(ID1) %>%
  get_summary_stats(height, type = "mean_sd")
ggboxplot(sp, x = "ID1", y = "height")
sp <- sp %>% mutate(TrialProv=factor(Trial):factor(Provenance))
sp <- sp %>% mutate(TrialProvFam=factor(Trial):factor(Provenance):factor(Family))
sp.trialProv <- sp %>%
  group_by(TrialProv) %>%
  get_summary_stats(height, type = "mean_sd")
ggboxplot(sp, x = "TrialProv", y = "height")
ggboxplot(sp, x = "TrialProvFam", y = "height")

# outliers
outliers <- sp %>% 
  group_by(TrialProvFam) %>%
  identify_outliers(height)
summary(factor(outliers$is.extreme)) # no extreme outliers

# check normality by examining residuals
# linear model
lm1  <- lm(height ~ TrialProvFam, data = sp)
# Create a QQ plot of residuals
ggqqplot(residuals(lm1))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(lm1))
# in the QQ plot, if the points fall approximately along the reference line, we can assume normality
# a little off at the extremes is ok but this seems more off...
# in the Shapiro-Wilk test, if the p-value is not significant, we can assume normality
# data is not normal...
hist(sp$height)
# if not normal, use kruskal-wallis
res.kruskalTPF <- sp %>% kruskal_test(height ~ TrialProvFam)
res.kruskalTPF
sp %>% kruskal_effsize(height ~ TrialProvFam)
# large effect size of trial:provenance:family (0.24)

library(corrplot)
# corrplot for provenance climate
corrplot(cor(sp[,c(11,16:35)]), method = "ellipse")
# corrplot for trial climate
corrplot(cor(sp[,c(11,36:55)]), method = "ellipse")

# based on Richards thesis (annual precipitation and growing season length important for clustering)
# fit mixed-effects model

mod1 <- lmer(height ~ MAP_T + (1|TrialProvFam), data = sp)
tidy(mod1)
# save model predictions as column to original data
sp$mod1<- predict(mod1)
summary(mod1)
# predicts increase in height with MAP_T but not significant

#include interaction
mod2 <- lmer(height ~ MAP_T + MAP_T*MAP_P + (1|TrialProvFam), data = sp)
tidy(mod2)
sp$mod2<- predict(mod2)
summary(mod2)
# predicts an increase in height with MAP_T and MAP_P but neither significant
# interaction predicts a decrease but not significant

mod3 <- lmer(height ~ MWMT_T + MWMT_T*MWMT_P + (1|TrialProvFam), data = sp)
summary(mod3)
# MWMT_T increases height, but not significant
# TrialProvFam explains 4% of the variance... (383/(383+7828))

mod4 <- lmer(height ~ MWMT_T + MWMT_T*MWMT_P +
               + PAS_T + PAS_T*PAS_P + (1|TrialProvFam), data = sp)
summary(mod4)
# simpsons paradox... predicted decrease in height with increasing MWMT_T when PAS is included
# predicts an increase in height with PAS_T
# no variable signficant again

mod5 <- lmer(height ~ CMD_T + Eref_T + TD_T + CMD_T*CMD_P + Eref_T*Eref_P + TD_T*TD_P + 
               (1|TrialProvFam), data = sp)
summary(mod5)
r2(mod5)
icc(mod5)
check_model(mod5)

compare_performance(mod1,mod2,mod3,mod4,mod5)
