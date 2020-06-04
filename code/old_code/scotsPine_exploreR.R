
library(tidyverse)
library(ggpubr) # publication ready plots
library(rstatix) # provides pipe-friendly R functions for easy statistical analyses

# read in raw data
wd <- "C:/Users/vanessa.burton/OneDrive - Forest Research/Documents/R/DeltaTraitSDM/" # PC wd
sp.raw <- read.csv(paste0(wd,"Scots_pine/Scots_pine_H.csv")) # raw data PC
sp.raw$X <- NULL
#str(sp.raw)
sp.raw<-na.omit(sp.raw)
summary(sp.raw)
colnames(sp.raw)[11]<-"height"

# summary stats
sp.prov <- sp.raw %>%
  group_by(ID1) %>%
  get_summary_stats(height, type = "mean_sd")
ggboxplot(sp.raw, x = "ID1", y = "height")

sp.trial <- sp.raw %>% 
  group_by(Trial) %>% 
  get_summary_stats(height, type = "mean_sd")
ggboxplot(sp.raw, x="Trial", y="height")

# both
sp.raw <- sp.raw %>% mutate(TrialProv=Trial:Provenance)
sp.tp <- sp.raw %>%
  group_by(TrialProv) %>%
  get_summary_stats(height, type = "mean_sd")
ggboxplot(sp.raw, x = "TrialProv", y = "height")

# outliers
outliers <- sp.raw %>% 
  group_by(ID1) %>%
  identify_outliers(height) # no extreme outliers

out_ids <- unique(outliers$id)

sp.nout <- sp.raw[ ! sp.raw$id %in% out_ids, ]
sp.tp2 <- sp.nout %>%
  group_by(TrialProv) %>%
  get_summary_stats(height, type = "mean_sd")
ggboxplot(sp.nout, x = "TrialProv", y = "height")

# check normality by examining residuals
# Build the linear model
model  <- lm(height ~ ID1, data = sp.nout)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
# In the QQ plot, if the points fall approximately along the reference line, we can assume normality.
# In the Shapiro-Wilk test, if the p-value is not significant, we can assume normality.
# my data is not normal...

# check normality by group
shapiro <- sp.nout %>%
  group_by(ID1) %>%
  shapiro_test(height)
ggqqplot(sp.nout, "height", facet.by = "ID1")

# if not normal, use kruskal-wallis
# plot residuals vs fit
#plot(model, 1)
# if no obvious relationship, then ok

# compute ANOVA
#res.aovP <- sp.nout %>% anova_test(height ~ ID1)
# the column ges corresponds to the generalized eta squared (effect size). 
# It measures the proportion of the variability in the outcome variable (here plant height) 
# that can be explained in terms of the predictor (here, provenance).
# An effect size of 0.015 (1.5%) means that 1.5% of the change in the height can be accounted for by provenance

# repeat for trial
#res.aovT <- sp.nout %>% anova_test(height ~ Trial)
# 12% change can be accounted for by trial

# compute Kruskal-Wallis
res.kruskalP <- sp.nout %>% kruskal_test(height ~ ID1)
res.kruskalP
sp.nout %>% kruskal_effsize(height ~ ID1)
# no significant differences between provenances with outliers, but when outliers are removed, there is a small significant effect (0.013)
res.kruskalT <- sp.nout %>% kruskal_test(height ~ Trial)
res.kruskalT
sp.nout %>% kruskal_effsize(height ~ Trial)
# large effect size of trial (0.21)

# A significant Kruskal-Wallis test is generally followed up by Dunn’s test to identify which groups are different.
# Pairwise comparisons
pwcT <- sp.nout %>% 
  dunn_test(height ~ Trial, p.adjust.method = "bonferroni") 
pwcT
# there is a statistically significant differences between trials 
# as assessed using the Kruskal-Wallis test (p = <0.0001). 
# Pairwise Wilcoxon test between groups showed that difference between trials were significant (p = <0.0001)

# visualise: box plots with p-values
pwcT <- pwcT %>% add_xy_position(x = "Trial")
ggboxplot(sp.nout, x = "Trial", y = "height") +
  stat_pvalue_manual(pwcT, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskalT, detailed = TRUE),
    caption = get_pwc_label(pwcT)
  )

# data not normal so don't use ANOVA
# two-way ANOVA
#plasticity <- sp.nout %>%
  #group_by(Trial, Provenance) %>%
  #get_summary_stats(height, type = "mean_sd")
#plasticity
# visualise
bxp <- ggboxplot(
  sp.nout, x = "Trial", y = "height",
  color = "Provenance")
bxp

# outliers
#outliers2 <- sp.nout %>%
  #group_by(Provenance, Trial) %>%
  #identify_outliers(height)
# some extreme outliers
#outliers2 %>% filter(is.extreme==TRUE)

# check normality
# Build the linear model with interactions between Trial and Provenance
model_int  <- lm(height ~ Trial*Provenance,
             data = sp.nout)
# Create a QQ plot of residuals
ggqqplot(residuals(model_int)) # not normal
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model_int))
ggqqplot(sp.nout, "height", ggtheme = theme_bw()) +
  facet_grid(Trial ~ Provenance)
# not normal
# so don't use anova
#res.aovTP <- sp.nout %>% anova_test(height ~ Provenance * Trial)
#res.aovTP
# There is a statistically significant interaction between Provenance and Trial for height,
# F(40, 1729) = 1.791, p = <0.05.
# A significant two-way interaction indicates that the impact that one factor (e.g., Provenance) 
# has on the outcome variable (e.g., height) depends on the level of the other factor (e.g., Trial)
# (and vice versa). 

# pairwise comparisons (between provenances)
library(emmeans)
pwcP <- sp.nout %>% 
  group_by(Trial) %>%
  emmeans_test(height ~ Provenance, p.adjust.method = "bonferroni") 
pwcP
# some significant differences between provenances, but mostly non-significant

emeans <- sp.nout %>% 
  emmeans_test(
    height ~ Provenance, p.adjust.method = "bonferroni",
    model = model_int
  )
summary(emeans)
emeans$p.adj.signif # most non-significant but some significant differences

# kruskal-wallis2
# compute Kruskal-Wallis
res.kruskalTP <- sp.nout %>% kruskal_test(height ~ TrialProv)
res.kruskalTP
sp.nout %>% kruskal_effsize(height ~ TrialProv)
# significant differences p<0.0001
# large effect size of trial:provenance (0.25)

# A significant Kruskal-Wallis test is generally followed up by Dunn’s test to identify which groups are different.
# Pairwise comparisons
pwcTP <- sp.nout %>% 
  dunn_test(height ~ TrialProv, p.adjust.method = "bonferroni") 
pwcTP
# there is a statistically significant differences between (some) provenances within trials 
# Pairwise Wilcoxon test between groups showed that difference between (some) provenances within trials were significant)

# visualise: box plots with p-values
pwcTP <- pwcTP %>% add_xy_position(x = "TrialProv")
ggboxplot(sp.nout, x = "TrialProv", y = "height") +
  coord_flip()+
  stat_pvalue_manual(pwcTP, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskalTP, detailed = TRUE),
    caption = get_pwc_label(pwcTP)
  )
