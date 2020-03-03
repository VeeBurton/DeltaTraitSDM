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
# family is a code for the provenance (mother trees known) and tree number. fathers not known
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
colnames(normal)[2]<-'Provenance'
colnames(sp)[8]<-'Provenance'
sp2<-merge(sp,normal[-c(22:24),], by='Provenance')
head(sp2)
head(sp2[,c(16:35)])
colnames(sp2)[16:35]<-c("MAT_P","MWMT_P","MCMT_P","TD_P","MAP_P","MSP_P","AHM_P","SHM_P","DD0_P","DD5_P","DD_18_P", "DD18_P", "NFFD_P","bFFP_P","eFFP_P","FFP_P","PAS_P",
                        "EMNT_P","Eref_P", "CMD_P")

normal$ID1<-as.character(normal$ID1)
normal$ID1[22]<-"BORDERS"
normal$ID1[23]<-"INVEREWE"
normal$ID1[24]<-"GLENSAUGH"
colnames(normal)[1]<-"Trial"

colnames(sp2)[4]<-"Trial"
sp3<-merge(sp2,normal[-c(1:21),], by='Trial')
head(sp3)
sp3<-sp3[,-c(36:39)]
head(sp3[,c(36:55)])
colnames(sp3)[36:55]<-c("MAT_T","MWMT_T","MCMT_T","TD_T","MAP_T","MSP_T","AHM_T","SHM_T","DD0_T","DD5_T","DD_18_T", "DD18_T", "NFFD_T","bFFP_T","eFFP_T","FFP_T","PAS_T",
                        "EMNT_T","Eref_T", "CMD_T")

summary(sp3)
colnames(sp3)[2]<-"Provenance"
colnames(sp3)[13:15]<-c("Latitude","Longitude","Elevation")
head(sp3)
#write.csv(sp3, "./Scots_pine/Scots_pine_H.csv")
sp3<-read.csv("./Scots_pine/Scots_pine_H.csv")
sp3$X<-NULL

# remove NAs
sp3<-na.omit(sp3)

# general exploratory plots
ggplot(sp3, aes(FFP_P,W17Height, colour=Provenance))+
  geom_point()+
  facet_wrap(~Trial)

par(mfrow=c(4,2),mar=c(3,3,3,1))
dotchart(sp3$W17Height, main="Height", group=sp3$Trial)
dotchart(sp3$Elevation, main="Elevation", group=sp3$Trial)
dotchart(sp3$MAT_P, main="MAT_P", group=sp3$Trial)
dotchart(sp3$MAT_T, main="MAT_T", group=sp3$Trial)
dotchart(sp3$DD5_P, main="DD5_P", group=sp3$Trial)
dotchart(sp3$DD5_T, main="DD5_T", group=sp3$Trial)
dotchart(sp3$FFP_P, main="FFP_P", group=sp3$Trial)
dotchart(sp3$FFP_T, main="FFP_T", group=sp3$Trial)
dev.off()

summary(sp3)

# before transforming
sp3$DD18_T<-NULL
corrplot(cor(sp3[,c(11,15:54)]), method = "ellipse")

# transform variables with large values as they will dominate any correlation coefficient
sp3$W17Height<-log(sp3$W17Height)
sp3$Elevation<-log(sp3$Elevation)
sp3$MAP_P<-log(sp3$MAP_P)
sp3$MSP_P<-log(sp3$MSP_P)
sp3$DD0_P<-log(sp3$DD0_P)
sp3$DD5_P<-log(sp3$DD5_P)
sp3$DD_18_P<-log(sp3$DD_18_P)
sp3$NFFD_P<-log(sp3$NFFD_P)
sp3$bFFP_P<-log(sp3$bFFP_P)
sp3$eFFP_P<-log(sp3$eFFP_P)
sp3$FFP_P<-log(sp3$FFP_P)
sp3$PAS_P<-log(sp3$PAS_P)
sp3$Eref_P<-log(sp3$Eref_P)
sp3$MAP_T<-log(sp3$MAP_T)
sp3$MSP_T<-log(sp3$MSP_T)
sp3$DD0_T<-log(sp3$DD0_T)
sp3$DD5_T<-log(sp3$DD5_T)
sp3$DD_18_T<-log(sp3$DD_18_T)
sp3$NFFD_T<-log(sp3$NFFD_T)
sp3$bFFP_T<-log(sp3$bFFP_T)
sp3$eFFP_T<-log(sp3$eFFP_T)
sp3$FFP_T<-log(sp3$FFP_T)
sp3$PAS_T<-log(sp3$PAS_T)
sp3$Eref_T<-log(sp3$Eref_T)

head(sp3[,c(11,15:54)])

# after transformation
corrplot(cor(sp3[,c(11,15:54)]), method = "ellipse") # no diff
pairs(sp3[,c(11,15:54)], na.action(na.omit))
collinearity <- as.data.frame(cor(sp3[,c(11,15:54)]))

cor <- collinearity %>% 
  mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>% 
  gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE) 

high.cor <- filter(cor, value <=-0.6 | value >=0.6)
ok.cor <- filter(cor, value >-0.6 & value <=0.6)

vars<-unique(ok.cor$Var2)
vars<-as.character(vars)

sp4 <- sp3[, (names(sp3) %in% vars)]

# Choosing variables
# Stepwise modelling
# Combination of variables (PCA etc)
# Look at variable colinearity

# look at distribution of height (response variable)
hist(sp3$W17Height) 

# standardise all explanatory variables (everything except Height and random effects)
#library(robustHD)
#head(sp3)
#head(sp3[,c(6:45)])
#sp3[,c(6:45)]<-robustHD::standardize(sp3[,c(6:45)], centerFun = mean)
#head(sp3)
#sp3$DD18_T<-NULL
#summary(sp3)
#sp3<-na.omit(sp3)

#sp3$sample<-as.character(sp3$sample)

# detect multicollinearity using VIF 
# split the data into training and test set
training.samples <- sp4$W17Height %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- sp4[training.samples, ]
test.data <- sp4[-training.samples, ]

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
sp4<-sp4[,-which(names(sp4) %in% ld.vars)] 

# variance and covariance
#####################################################

str(sp4)
#sp4.s<-sp4[,-c(3,4,7,50)] # remove factors
var(sp4)
cor(sp4)
corrplot(cor(sp4), method = "ellipse")

training.samples <- sp4$W17Height %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- sp4[training.samples, ]
test.data <- sp4[-training.samples, ]

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

# detect multicollinearity
VIF <- car::vif(model2) %>%
  as.list() %>% 
  as.data.frame() %>% 
  gather(key = 'variable', value = 'VIF') %>% 
  arrange(desc(VIF))

# remove variables with high VIF (above 5-10)
best.vars<-unique(VIF$variable[which(VIF$VIF<100)])
best.vars



####################################################################################################
# good practice to standardise all explanatory variables (x)
# include Trial/Block/Tree_id and Provenance as random effects - (1|Trial/Block/Tree_id)
# may need to make nesting explicit by creating new variables e.g. Trial1-Blocka, Trial1-Blockb etc.

# may need to make nesting explicit by creating new variables e.g. Trial1-Blocka, Trial1-Blockb etc.
# make nested variables
# site/block/population/family/seedling (or tag?)
head(sp3)
sp3 <- within(sp3, sample <- factor(Trial:Provenance:Family:Block))
sample <- unique(sp3$sample)
summary(sample)
length(unique(sp3$sample))

sp_summary<-sp3 %>% 
  group_by(Trial) %>% 
  summarise(Blocks = length(unique(Block)),
            Provenances = length(unique(Provenance)),
            Families = length(unique(Family)),
            Seedlings = length(unique(Seedling)),
            Individuals = length(unique(Tag)),
            Observations = sum(n()))

head(sp3[,-c(1:10)])
sp3<-sp3[,-c(1:10)]

# plot data by trial/provenance etc.
ggplot(aes(W17Height), data = sp3) + geom_histogram(binwidth = 40) +
  facet_wrap(~ Trial) +
  xlab("Height") + ylab("Frequency")

ggplot(aes(W17Height), data = sp3) + geom_histogram(binwidth = 40) +
  facet_wrap(~ Provenance) +
  xlab("Height") + ylab("Frequency")

boxplot(W17Height ~ Trial, data = sp3)
boxplot(W17Height ~ Provenance, data = sp3)

pinus <- sp4
pinus$Trial<- sp3$Trial
pinus$Block <- sp3$Block

# e.g.
SPmod1 <- lmer(W17Height ~ DD_18_T + MAT_T + (1|Trial/Block), data = pinus)

summary(SPmod1)
coef(SPmod1)
augment(SPmod1)

# plot residuals
plot(SPmod1, which = 1)
# q plot
qqnorm(resid(SPmod1))
qqline(resid(SPmod1))

# not a good model! need more vars and to work out why i'm getting so many aliased coefficients
# job for next week!

