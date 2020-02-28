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
library(broom)
library(tidyverse)
library(corrplot)
library(caret)
library(ggplot2)
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

# height in mm
# there are 3 sites – called Borders (Yair), Glensaugh, Inverewe here – and each is a fully randomised block design: 
# there are 4 blocks at Borders, Glensaugh and 3 at Inverewe
# in each site there are 21 PlantingSites, 8 families and either 3 or 4 individuals per family (1 per block)
# so there are 21*8*4=672 at Borders, Glensaugh; 21*8*3= 504 at Inverewe; total = 1848
# each tree has a unique identity (Tag), and details of each Seed Zone, PlantingSite, Family and Individual number are given per tree.
# ‘PlantingSite’ tallies with the sites of origin in the coordinate file I sent you previously.

#locations <- read.csv("./Scots_pine/scots_pine_source_locations.csv")
#head(locations)
#colnames(locations)[1]<-"Name"
#codes<-as.data.frame(unique(sp$Population))
#colnames(codes)[1]<-'code'
#codes<-codes[order(codes$code),]
#locations$Population<-codes
#locations<-data.frame(locations$Name, locations$Population, locations$X, locations$Y)
#colnames(locations)<-c('ID1', 'ID2', 'Lat', 'Long')
# need elev for Climate EU
#write.csv(locations, "./Scots_pine/scots_pine_source_locations2.csv")
# extracted elev from 50m DEM in Arc, for trial sites and provenances/sources - running through Climate EU
all.locations <- read.csv("./Scots_pine/scots_pine_all_locations_elev.csv")
ggplot(aes(Long,Lat, color=ID1), data=all.locations)+geom_point()

normal<-read.csv("./Scots_pine/scots_pine_all_locations_elev_Normal_1961_1990Y.csv")
s2050<-read.csv("./Scots_pine/scots_pine_all_locations_elev_HadGEM2-ES_rcp85_2050sY.csv")
s2080<-read.csv("./Scots_pine/scots_pine_all_locations_elev_HadGEM2-ES_rcp85_2080sY.csv")

# merge to trait data
colnames(normal)[2]<-'Population'
sp2<-merge(sp,normal[-c(22:24),], by='Population')
head(sp2)
head(sp2[,c(16:35)])
colnames(sp2)[16:35]<-c("MAT_P","MWMT_P","MCMT_P","TD_P","MAP_P","MSP_P","AHM_P","SHM_P","DD0_P","DD5_P","DD_18_P", "DD18_P", "NFFD_P","bFFP_P","eFFP_P","FFP_P","PAS_P",
                        "EMNT_P","Eref_P", "CMD_P")

normal$ID1<-as.character(normal$ID1)
normal$ID1[22]<-"BORDERS"
normal$ID1[23]<-"INVEREWE"
normal$ID1[24]<-"GLENSAUGH"
colnames(normal)[1]<-"PlantingSite"

sp3<-merge(sp2,normal[-c(1:21),], by='PlantingSite')
head(sp3)
sp3<-sp3[,-c(36:39)]
head(sp3[,c(36:55)])
colnames(sp3)[36:55]<-c("MAT_T","MWMT_T","MCMT_T","TD_T","MAP_T","MSP_T","AHM_T","SHM_T","DD0_T","DD5_T","DD_18_T", "DD18_T", "NFFD_T","bFFP_T","eFFP_T","FFP_T","PAS_T",
                        "EMNT_T","Eref_T", "CMD_T")

summary(sp3)
colnames(sp3)[2]<-"Population"
colnames(sp3)[13:15]<-c("Latitude","Longitude","Elevation")
head(sp3)
#write.csv(sp3, "./Scots_pine/Scots_pine_H.csv")
sp3<-read.csv("./Scots_pine/Scots_pine_H.csv")
sp3$X<-NULL

# remove NAs
sp3<-na.omit(sp3)

# may need to make nesting explicit by creating new variables e.g. Trial1-Blocka, Trial1-Blockb etc.
# make nested variables
# site/block/population/family/seedling (or tag?)
sp3 <- within(sp3, sample <- factor(PlantingSite:Block:Population))
head(sp3)
sp3<-sp3[,-c(1,2,5,6,7)]

# plot data by trial/provenance etc.
ggplot(aes(W17Height), data = sp3) + geom_histogram(binwidth = 40) +
  facet_wrap(~ PlantingSite) +
  xlab("Height") + ylab("Frequency")

ggplot(aes(W17Height), data = sp3) + geom_histogram(binwidth = 40) +
  facet_wrap(~ Population) +
  xlab("Height") + ylab("Frequency")

boxplot(W17Height ~ PlantingSite, data = sp3)
boxplot(W17Height ~ Population, data = sp3)

# Choosing variables
# Stepwise modelling
# Combination of variables (PCA etc)
# Look at variable colinearity

# look at distribution of height (response variable)
hist(sp3$W17Height) 

# standardise all explanatory variables (everything except Height and random effects)
library(robustHD)
head(sp3)
head(sp3[,c(8:50)])
sp3[,c(8:50)]<-robustHD::standardize(sp3[,c(8:50)], centerFun = mean)
head(sp3)
sp3$DD18_T<-NULL
summary(sp3)
sp3<-na.omit(sp3)

# detect multicollinearity using VIF 
# split the data into training and test set
training.samples <- sp3$W17Height %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- sp3[training.samples, ]
test.data <- sp3[-training.samples, ]

# build a regression model with all variables
model1 <- lm(W17Height ~., data = train.data)
# make predictions
predictions <- model1 %>% predict(test.data)
# model performance
data.frame(
  RMSE = RMSE(predictions, test.data$W17Height),
  R2 = R2(predictions, test.data$W17Height)
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
ldvars<-ld.vars[48:88]
sp3<-sp3[,-which(names(sp3) %in% ldvars)] # this is all explanatory variables :(

# variance and covariance
#####################################################

str(sp3)
sp3$id<-as.numeric(sp3$id)
sp3$Tag<-as.numeric(sp3$Tag)
sp3$Seedling<-as.numeric(sp3$Seedling)
sp3$W17Height<-as.numeric(sp3$W17Height)
str(sp3)

sp3.s<-sp3[,-c(3,4,7,50)] # remove factors
var(sp3.s)
cor(sp3.s)
corrplot(cor(sp3.s), method = "ellipse")

training.samples <- sp3.s$W17Height %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- sp3.s[training.samples, ]
test.data <- sp3.s[-training.samples, ]

# build a regression model with all variables
model2 <- lm(W17Height ~., data = train.data)
# make predictions
predictions <- model2 %>% predict(test.data)
# model performance
data.frame(
  RMSE = RMSE(predictions, test.data$W17Height),
  R2 = R2(predictions, test.data$W17Height)
)
model2_summary<-tidy(model2)

# detect multicollinearity
car::vif(model2) 

# just trial vars
head(sp3)
sp_T<-sp3[,c(50,1:10,31:49)]
str(sp_T)
sp_T2<-sp_T[,-c(1,4,5,8)] # remove factors

corrplot(cor(sp_T2), method = "ellipse")
cor_spT<-as.data.frame(cor(sp_T2))
cor_spT[which (cor_spT < 0.4),]

sp_T$Latitude<-NULL
sp_T$Longitude<-NULL

training.samples <- sp_T$W17Height %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- sp_T[training.samples, ]
test.data <- sp_T[-training.samples, ]

# build a regression model with all variables
modelT <- lm(W17Height ~., data = train.data)
# make predictions
predictions <- modelT %>% predict(test.data)
# model performance
data.frame(
  RMSE = RMSE(predictions, test.data$W17Height),
  R2 = R2(predictions, test.data$W17Height)
)
modelT_summary<-tidy(modelT)

# detect multicollinearity
car::vif(modelT) 
summary(modelT)$coeff


ld.vars <- attributes(alias(model2)$Complete)$dimnames[[1]]
ld.vars
sp3<-sp3[,-which(names(sp3) %in% ldvars)]

# detect multicollinearity
VIF <- car::vif(model2) %>%
  as.list() %>% 
  as.data.frame() %>% 
  gather(key = 'variable', value = 'VIF') %>% 
  arrange(desc(VIF))

# remove variables with high VIF (above 5-10)
best.vars<-unique(VIF$variable[which(VIF$VIF<12)])
best.vars


####################################################################################################
# good practice to standardise all explanatory variables (x)
# include Trial/Block/Tree_id and Provenance as random effects - (1|Trial/Block/Tree_id)
# may need to make nesting explicit by creating new variables e.g. Trial1-Blocka, Trial1-Blockb etc.

# e.g.
SPmod1 <- lmer(W17Height ~ FFP_T
               + (1|PlantingSite/Block/Population), data = sp3)

summary(SPmod1)
coef(SPmod1)
augment(SPmod1)

# plot residuals
plot(SPmod1, which = 1)
# q plot
qqnorm(resid(SPmod1))
qqline(resid(SPmod1))

