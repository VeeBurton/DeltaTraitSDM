library(lme4)
library(nlme)
library(broom)
library(ggplot2)
library(performance)

### Multiple and Logistic Regression

Pinus<-read.csv("./Scots_pine/Scots_pine_H.csv")
Pinus$X<-NULL
# remove NAs
Pinus<-na.omit(Pinus)
colnames(Pinus)[11]<-'H'
head(Pinus)

Pinus <- Pinus %>% mutate(block=factor(Trial:Block),
                    popSite=factor(Trial:Provenance),
                    famSite=factor(Trial:Family))

# notes from datacamp course
# can include variables as both fixed and random effects
# if i've understood this correctly...
# fixed effect slope would estimate how height is changing across trials
# random effect slope corrects for trials having different changes in height

# can compare models with anova
# test if one model explains more variability than the other model
# If you wanted to see if trial site is important for predicting H , we can build a null model with only trial as a random-effect 
# and a trial model that includes Trial You can then compare the two models using the anova() function
# could do same for provenance

# Build the Null model with only provenance as a random-effect
null_model <- lmer(log(H) ~ (1 | Provenance) , data = Pinus)

# Build the Trial model with Trial as a fixed and random slope and provenance as the random-effect
trial_model <- lmer(log(H) ~ Trial + (1 + Trial | Provenance) , data = Pinus)
block_model <- lmer(log(H) ~ Block + (1 + Block | popSite), data = Pinus) # use explicitly nested vars (Block = Trial:Block, popSite=Trial:Provenance)

# Compare null_model and year_model using an anova
anova(null_model, trial_model)
anova(trial_model,block_model)

#############
# interaction terms
#############

# to include an interaction term
# the colon : represents the interaction
interaction_mod <- lm(H ~ factor(DD_18_T) + Trial + DD_18_T:Trial ,data = Pinus)
interaction_mod

ggplot(Pinus, aes(DD_18_P,H,colour=Trial))+
  geom_point()+
  geom_smooth(method='lm',se=FALSE)

ggplot(Pinus, aes(MAT_P,H,colour=Trial))+
  geom_point()+
  geom_smooth(method='lm',se=FALSE)

#############
# adding another numerical explanatory variable
#############
# data space is now 3D so ggplot won't work (no z aes)

# can tile the plane
grid<-Pinus %>% modelr::data_grid(mwmt=modelr::seq_range(MWMT_T,by=1),
                    pas=modelr::seq_range(PAS_T,by=1))
mod <- lm(H ~ MWMT_T + PAS_T, data=Pinus)
H_hats <- augment(mod,newdata = grid)

# or 3D visualisation
library(plotly)
x<-unique(Pinus$PAS_T)
y<-unique(Pinus$H)
plane
plot_ly(data=Pinus, z=~MWMT_T, x=~PAS_T, y=~H, opacity=0.6) %>% 
  add_markers(marker=list(size=2)) %>% 
  add_surface(x=~x,y=~y,z=~z, showcale=FALSE,cmax=1)

##############
# 26/03/2020
##############

# bin the data
ggplot(data = Pinus, mapping = aes(x=DD_18_T)) + 
  geom_histogram(aes(y=..density..),fill="bisque",color="white",alpha=0.7) + 
  geom_density() +
  geom_rug() +
  labs(x='Degree days above 18 for trial sites') +
  theme_minimal()
summary(Pinus$DD_18_T)
breaks <- c(3600,3650,3700,3750,3800,3850,3900,3950)
tags <- c('[3600-3650)','[3650-3700)','[3700-3750)','[3750-3800)','[3800-3850)','[3850-3900)','[3900-3950)')
group_tags <- cut(Pinus$DD_18_T, 
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)

# inspect bins
summary(group_tags)
# ordered factor
DD18T_groups <- factor(group_tags,
                       levels = tags,
                       ordered = TRUE)
ggplot(data = as_tibble(group_tags), mapping = aes(x=value)) + 
  geom_bar(fill="bisque",color="white",alpha=0.7) + 
  stat_count(geom="text", aes(label=sprintf("%.4f",..count../length(group_tags))), vjust=-0.5) +
  labs(x='DD18 for trial sites') +
  theme_minimal() 
Pinus$DD18bin <- factor(group_tags,
                        levels = tags,
                        ordered = TRUE)
data_space <- ggplot(data=Pinus,aes(x=DD18bin,y=W17Height))+
  geom_point()+
  geom_line()
data_space


# model (probability approach)
trial_lm <- lm(W17Height~DD18bin+Trial, data=Pinus)
summary(trial_lm)
exp(coef(trial_lm))
# augmented model
augie <- trial_lm %>% augment(type.predict='response')
# logistic model (probability)
data_space +
  geom_line(data=augie, aes(x=DD18bin,y=.fitted),color='red')


# odds approach
# compute odds for bins
Pinus <- Pinus %>% mutate(odds=W17Height/(1-W17Height))
# plot binned odds
data_space <- ggplot(data=Pinus,aes(DD18bin,odds))+geom_point()+geom_line()
# compute odds for observations
augie2 <- augie %>% mutate(odds_hat=.fitted/(1-.fitted))
# logistic model on odds scale
data_space +
  geom_line(data=augie2,aes(DD18bin,odds_hat), color = "red")

