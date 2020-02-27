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

sp <- read.csv("./Scots_pine/PS_Common_Garden_Heights.csv")
head(sp)
summary(sp)
str(sp)
colnames(sp)[1]<-"id"
summary(sp$W17Height)
sp$W17Height[which(sp$W17Height=="*")]<-NA
sp$W17Height<-as.numeric(sp$W17Height)
sp$Seedling[which(sp$Seedling=="*")]<-NA
sp$Seedling<-as.numeric(sp$Seedling)
summary(sp$Block)
sp$Block<-as.character(sp$Block)
sp$Block[which(sp$Block=="c")]<-"C"
sp$Block<-as.factor(sp$Block)

# there are 3 sites – called Borders, Glensaugh, Inverewe here – and each is a fully randomised block design: 
# there are 4 blocks at Borders, Glensaugh and 3 at Inverewe
# in each site there are 21 PlantingSites, 8 families and either 3 or 4 individuals per family (1 per block)
# so there are 21*8*4=672 at Borders, Glensaugh; 21*8*3= 504 at Inverewe; total = 1848
# each tree has a unique identity (Tag), and details of each Seed Zone, PlantingSite, Family and Individual number are given per tree.
# ‘PlantingSite’ tallies with the sites of origin in the coordinate file I sent you previously.

locations <- read.csv("./Scots_pine/scots_pine_source_locations.csv")
head(locations)
colnames(locations)[1]<-"Name"
codes<-as.data.frame(unique(sp$Population))
colnames(codes)[1]<-'code'
codes<-codes[order(codes$code),]
locations$Population<-codes
locations<-data.frame(locations$Name, locations$Population, locations$X, locations$Y)
colnames(locations)<-c('ID1', 'ID2', 'Lat', 'Long')
# need elev for Climate EU
write.csv(locations, "./Scots_pine/scots_pine_source_locations2.csv")

# plot data by trial/provenance etc.
ggplot(aes(W17Height), data = sp) + geom_histogram() +
  facet_wrap(~ PlantingSite) +
  xlab("Height") + ylab("Frequency")

boxplot(W17Height ~ PlantingSite, data = sp)

# Choosing variables
# Stepwise modelling
# Combination of variables (PCA etc)
# Look at variable colinearity

# look at distribution of height (response variable)
hist(sp$W17Height) 

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

