
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)

sp <- read.csv("~/R/DeltaTraitSDM/Scots_pine/sp_data_height_climateTP.csv")
head(sp)
sp<-na.omit(sp)

# order trials by longitude
sp$Trial <- ordered(sp$Trial, levels = c("INVEREWE", "BORDERS", "GLENSAUGH"))
# provenance boxplot by longitude
ggplot(sp)+
  geom_boxplot(aes(Longitude,height,color=ID1))+
  facet_wrap(~Trial)

############################################################################################################
# Perry et al. 2016
# To test for significant differences in height among families, populations and blocks
# nested analysis of variance (ANOVA) tests were performed with population as a fixed effect, 
# and families nested within population and block as random effects 
############################################################################################################

# test for significant differences in height among populations and trials
# nested anova, with Provenance as a fixed effect and Provenance nested within Trial as random effects
mod1 <- lmer(height ~ Provenance + (1|Trial/Provenance), sp, REML = TRUE)
summary(mod1)
hist(residuals(mod1), col="darkgray") # residuals are normally distributed
plot(fitted(mod1),  residuals(mod1))

# % variance
tot.var <- 474.2+105190+79261.4
trial.var <- 105190/tot.var*100 #56%
prov.var <- 474.2/tot.var*100 #0.25%
resid.var <- 79261.4/tot.var*100 #42%

# nested analysis of variance
# are the variations between provenance means due to true differences about the populations means or just due to sampling variability
# F statistics = Variation among sample means / Variation within groups
aov1 <- anova(mod1) # significant p <0.001 but small effect size F 7.2
# variation of height means between different provenances is slightly larger than the variation of height within each provenance,
rand(mod1)
difflsmeans(mod1, test.effs="Provenance") # least square means

# posthoc tests see *which* provenances are significantly different
library(multcomp)
posthoc <- glht(mod1,linfct = mcp(Provenance="Tukey"))
mcs <- summary(posthoc, test=adjusted("single-step"))
mcs
# plot
plot(posthoc)

############################################################################################################
# parallel slopes (DataCamp Multiple and Logistic Regression)
# understand how height varies as a function of climate/transfer distance (for temperature and precipitation) but also trial
############################################################################################################
library(broom)

sp2 <- sp %>% 
  mutate(CD_temp=MAT_T-MAT_P,
         CD_prec=MAP_T-MAP_P)
# CD = climate distance
# positive difference, tree moved to warmer/wetter temps
# negative difference, tree moved to colder/drier temps

# parallel slopes model
lm1 <- lm(height ~ CD_temp + factor(Trial), sp2)
summary(lm1)
lm1aug <- augment(lm1)
data_space <- ggplot(lm1aug, aes(x=CD_temp,y=height, color=factor.Trial.))+
  geom_point()
data_space +
  geom_line(aes(y=.fitted))+
  theme_minimal()
# interpretation
levels(sp2$Trial) # by default the reference (intercept) is the first level alphabetically - i can set this to change
# intercepts are in units of respnse, slope is in unit of response per unit of explanatory variable
# intercept : in units of response (height) - expected height for first trial (Inverewe?), CD of 0
# factor(Trial).L - this trial has height typically 59mm more than the first trial (intercept), after controlling for CD
# factor(Tral).Q - this trial has height typically 390mm less than the first trial (intercept), after control for CD
# CD estimate is our slope coefficient, each degree change (increase?) in climate distance is associated with a decrease in height of 76mm
# so. provenances moved to cooler climates grow taller?

# interaction model - allow different slopes
# add interaction manually
lm1i <- lm(height ~ CD_temp + Trial + CD_temp:Trial, sp2)
summary(lm1i)
# interpretation becomes more complicated
lm1i_aug <- augment(lm1i)
data_space <- ggplot(lm1i_aug, aes(x=CD_temp,y=height, color=Trial))+
  geom_point()
data_space +
  geom_line(aes(y=.fitted))+
  theme_minimal()
# drop in height with unit increase in climate distance
# larger rate of decrease in height the further west the seeds are transferred?
# CD:Inverewe -74 (drop in height per unit CD) -- most westerly trial
# CD:Borders -38 (drop in height per unit CD) -- middle trial
# CD:Glensuagh -10 (drop in height per unit CD) - most eastern trial
# this assumes the reference is Inverewe and Trial.L is borders and Trial.Q Glensaugh

# provenance interaction model
provInt <- lm(height ~ CD_temp + Provenance + CD_temp:Provenance, sp2)
summary(provInt)
provIntAug <- augment(provInt)
data_space <- ggplot(provIntAug, aes(x=CD_temp,y=height, color=Provenance))+
  geom_point()
data_space +
  geom_line(aes(y=.fitted))+
  theme_minimal()+
  xlab("Transfer Distance (MAT)")

ggplot(data= provIntAug,aes(x=.fitted,y=height))+
  geom_point()+
  geom_smooth(method = "lm", se=0)

# repeat for precipitation
lm2 <- lm(height ~ CD_prec + factor(Trial), sp2)
summary(lm2)
lm2aug <- augment(lm2)

data_space <- ggplot(lm2aug, aes(x=CD_prec,y=height, color=factor.Trial.))+
  geom_point()
data_space +
  geom_line(aes(y=.fitted))+
  theme_minimal()+
  xlab("Transfer distance (MAP)")

# considers size of transfer distance while also accounting for different Trial conditions

############################################################################################################
# George et al. 2020
# fitted a general linear model between height and distance between trial site climate and provenance climate
# transfer distance = trial climate - prov climate
############################################################################################################

# 1. assigned provenances to climatic clusters based on PCA
 
# 2. test for general patterns of local adaptation
# general linear model
# separately for each trial site
# height ~ MAT_diff + MAP_diff

borders <- sp %>% 
  filter(Trial=="BORDERS") %>% 
  mutate(MAT_diff=MAT_T-MAT_P,
         MAP_diff=MAP_T-MAP_P)
# negative differences mean tree moved to colder/wetter conditions
# positive differences mean tree moved to warmer/drier conditions
glmB <- glm(height ~ MAT_diff + MAP_diff, borders, family=gaussian)
summary(glmB)

inverewe <- sp %>% 
  filter(Trial=="INVEREWE") %>% 
  mutate(MAT_diff=MAT_T-MAT_P,
         MAP_diff=MAP_T-MAP_P)
glmI <- glm(height ~ MAT_diff + MAP_diff, inverewe, family=gaussian)
summary(glmI)

glensaugh <- sp %>% 
  filter(Trial=="GLENSAUGH") %>% 
  mutate(MAT_diff=MAT_T-MAT_P,
         MAP_diff=MAP_T-MAP_P)
glmG <- glm(height ~ MAT_diff + MAP_diff, glensaugh, family=gaussian)
summary(glmG)

# 3. Test for variation in phenotypic plasticity of growth traits (GxE) and to estimate its contribution to overall phenotypic variation 
# formulated a mixed model
# Yijklm = β0 + β1Si + β2Pj + β3Mk(Pj) + β4Bl(Si) + β5(SiPj) + β6(SiMk(Pj)) + eijklm
# With Yijklm being the phenotype of the mth tree, belonging to the lth block nested within the ith site (Bl(Si)), 
# belonging to the kth mother tree nested within provenance j (Mk(Pj), originating from provenance j (Pj) and growing in the ith trial site (Si). 
# eijklm is a random error term and SiPj and SiMk(Pj) are the crossed genotype x environment interaction terms separated for provenance-by- site 
# and family-by-site interactions, respectively.
