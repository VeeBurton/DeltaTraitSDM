
library(tidyverse)
library(ggpubr) # publication ready plots
library(rstatix) # provides pipe-friendly R functions for easy statistical analyses

# read in raw data
sp.raw <- read.csv(paste0(wd,"Scots_pine/Scots_pine_H.csv")) # raw data PC
sp.raw$X <- NULL
#str(sp.raw)
sp.raw<-na.omit(sp.raw)
summary(sp.raw)
colnames(sp.raw)[11]<-"H"

# summary stats
sp.prov <- sp.raw %>%
  group_by(ID1) %>%
  get_summary_stats(H, type = "mean_sd")
ggboxplot(sp.raw, x = "ID1", y = "H")

sp.trial <- sp.raw %>% 
  group_by(Trial) %>% 
  get_summary_stats(H, type = "mean_sd")
ggboxplot(sp.raw, x="Trial", y="H")

# outliers
outliers <- sp.raw %>% 
  group_by(ID1) %>%
  identify_outliers(H) # no extreme outliers

# check normality by examining residuals
# Build the linear model
model  <- lm(H ~ ID1, data = sp.raw)
# Create a QQ plot of residuals
ggqqplot(residuals(model))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model))
# In the QQ plot, if the points fall approximately along the reference line, we can assume normality.
# In the Shapiro-Wilk test, if the p-value is not significant, we can assume normality.
# my data is not normal...

# check normality by group
shapiro <- sp.raw %>%
  group_by(ID1) %>%
  shapiro_test(H)
ggqqplot(sp.raw, "H", facet.by = "ID1")

# if not normal, use kruskal-wallis

# plot residuals vs fit
plot(model, 1)
# if no obvious relationship, then ok

# compute ANOVA
res.aov1 <- sp.raw %>% anova_test(H ~ ID1)
# the column ges corresponds to the generalized eta squared (effect size). 
# It measures the proportion of the variability in the outcome variable (here plant height) 
# that can be explained in terms of the predictor (here, provenance).
# An effect size of 0.015 (1.5%) means that 1.5% of the change in the height can be accounted for by provenance

# repeat for trial
res.aov2 <- sp.raw %>% anova_test(H ~ Trial)
# 12% change can be accounted for by trial

# compute Kruskal-Wallis
res.kruskal1 <- sp.raw %>% kruskal_test(H ~ ID1)
res.kruskal1
sp.raw %>% kruskal_effsize(H ~ ID1)
# no significant differences between provenances
res.kruskal2 <- sp.raw %>% kruskal_test(H ~ Trial)
res.kruskal2
sp.raw %>% kruskal_effsize(H ~ Trial)
# large effect size of trial (0.19)

# A significant Kruskal-Wallis test is generally followed up by Dunn’s test to identify which groups are different.
# Pairwise comparisons
pwc <- sp.raw %>% 
  dunn_test(H ~ Trial, p.adjust.method = "bonferroni") 
pwc

# There was a statistically significant differences between trials 
# as assessed using the Kruskal-Wallis test (p = <0.0001). 
# Pairwise Wilcoxon test between groups showed that difference between trials were significant (p = <0.0001)

# visualise: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Trial")
ggboxplot(sp.raw, x = "Trial", y = "H") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal2, detailed = TRUE),
    caption = get_pwc_label(pwc)
  )

# two-way ANOVA
plasticity <- sp.raw %>%
  group_by(Trial, Provenance) %>%
  get_summary_stats(H, type = "mean_sd")
plasticity
# visualise
bxp <- ggboxplot(
  sp.raw, x = "Trial", y = "H",
  color = "Provenance")
bxp

# outliers
outliers2 <- sp.raw %>%
  group_by(Provenance, Trial) %>%
  identify_outliers(H)
# some extreme outliers
outliers2 %>% filter(is.extreme==TRUE)

# check normality
# Build the linear model
model2  <- lm(H ~ Provenance*Trial,
             data = sp.raw)
# Create a QQ plot of residuals
ggqqplot(residuals(model2))
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model2))
ggqqplot(sp.raw, "H", ggtheme = theme_bw()) +
  facet_grid(Provenance ~ Trial)
# not normal

# try anova anyway
res.aov3 <- sp.raw %>% anova_test(H ~ Provenance * Trial)
res.aov3
# There is a statistically significant interaction between Provenance and Trial for height,
# F(40, 1729) = 1.791, p = <0.05.
# A significant two-way interaction indicates that the impact that one factor (e.g., Provenance) 
# has on the outcome variable (e.g., height) depends on the level of the other factor (e.g., Trial)
# (and vice versa). 

# pairwise comparisons (between provenances)
library(emmeans)
pwc2 <- sp.raw %>% 
  group_by(Provenance) %>%
  emmeans_test(H ~ Trial, p.adjust.method = "bonferroni") 
pwc2
# some significant differences between provenances

emeans <- sp.raw %>% 
  emmeans_test(
    H ~ Trial, p.adjust.method = "bonferroni",
    model = model2
  )
summary(emeans)
emeans$p.adj.signif

pwc2 <- pwc2 %>% add_xy_position(x = "Trial")
bxp +
  stat_pvalue_manual(pwc2) +
  labs(
    subtitle = get_test_label(res.aov3, detailed = T),
    caption = get_pwc_label(pwc2)
  )

# kruskal-wallis2
# compute Kruskal-Wallis
# explicitly nest
sp.raw <- sp.raw %>% mutate(TP=Provenance:Trial)
res.kruskal2 <- sp.raw %>% kruskal_test(H ~ TP)
res.kruskal2
sp.raw %>% kruskal_effsize(H ~ TP)
# significant differences
# large effect size of trial:provenance (0.21)

# A significant Kruskal-Wallis test is generally followed up by Dunn’s test to identify which groups are different.
# Pairwise comparisons
pwc3 <- sp.raw %>% 
  dunn_test(H ~ TP, p.adjust.method = "bonferroni") 
pwc3

# There was a statistically significant differences between (some) provenances within trials 
# Pairwise Wilcoxon test between groups showed that difference between (some) provenances within trials were significant)

# visualise: box plots with p-values
pwc3 <- pwc3 %>% add_xy_position(x = "TP")
ggboxplot(sp.raw, x = "TP", y = "H") +
  stat_pvalue_manual(pwc3, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.kruskal2, detailed = TRUE),
    caption = get_pwc_label(pwc3)
  )
