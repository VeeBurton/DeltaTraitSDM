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
library(performance)
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

sp_summary<-sp3 %>% 
  group_by(Trial) %>% 
  summarise(Blocks = length(unique(Block)),
            Provenances = length(unique(Provenance)),
            Families = length(unique(Family)),
            Seedlings = length(unique(Seedling)),
            Individuals = length(unique(Tag)),
            Observations = sum(n()))

# general exploratory plots
ggplot(sp3, aes(FFP_P,W17Height, colour=Provenance))+
  geom_point()+
  facet_wrap(~Trial)


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

cor$combo<-paste0(cor$Var1,"-",cor$Var2)
cor<-cor[unique(cor$combo),]
head(cor)
cor<- cor[order(cor$value),]

high.cor <- filter(cor, value <=-0.6 | value >=0.6)
ok.cor <- filter(cor, value >-0.2 & value <=0.2)
ok.cor<- ok.cor[order(ok.cor$value),]

vars<-unique(ok.cor$combo)
vars<-as.character(vars)
vars

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
best.vars<-unique(VIF$variable[which(VIF$VIF<10)])
best.vars

####################################################################################################
# good practice to standardise all explanatory variables (x)
# include Trial/Block/Tree_id and Provenance as random effects - (1|Trial/Block/Tree_id)
# may need to make nesting explicit by creating new variables e.g. Trial1-Blocka, Trial1-Blockb etc.

# may need to make nesting explicit by creating new variables e.g. Trial1-Blocka, Trial1-Blockb etc.
# make nested variables
# site/block/population/family/seedling (or tag?)


# pinus <- read.csv("./Scots_pine/Scots_pine_H_standardised_allvars.csv")
# just variables identified by PCA (along same dimension as height)
PCA_vars<-c("W17Height","MAT_T","MAP_T","FFP_T","DD5_T","MSP_T","eFFP_T","NFFD_T","MCMT_T","EMNT_T","MWMT_T","AHM_T","SHM_T",
            "TD_T","DD0_T","CMD_T","DD_18_T","Eref_T","bFFP_T")
pinus<-pinus[,which(names(pinus) %in% PCA_vars)] 
pinus$Trial<- sp3$Trial
pinus$Block <- sp3$Block
pinus$Family <- sp3$Family
pinus$Provenance <- sp3$Provenance
pinus <- within(pinus, sample <- factor(Provenance:Family:Block))
head(pinus)
pinus<-pinus[,-c(21:23)]

boxplot(W17Height ~ Trial, data = pinus)
boxplot(W17Height ~ sample, data = pinus)

corrplot(cor(pinus[,-c(20:21)]), method = "ellipse")
pairs(pinus[,-c(20:21)])

training.samples <- pinus$W17Height %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- pinus[training.samples, ]
test.data <- pinus[-training.samples, ]

# build a regression model with all variables
model3 <- lm(W17Height ~., data = train.data)
# make predictions
predictions <- model3 %>% predict(test.data)
# model performance
data.frame(
  RMSE = RMSE(predictions, test.data$W17Height),
  R2 = R2(predictions, test.data$W17Height)
)
model3_summary<-tidy(model3)

ld.vars <- attributes(alias(model3)$Complete)$dimnames[[1]]
ld.vars
pinus2<-pinus[,-which(names(pinus) %in% ld.vars)] 

training.samples <- pinus2$W17Height %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- pinus2[training.samples, ]
test.data <- pinus2[-training.samples, ]

# build a regression model with all variables
model4 <- lm(W17Height ~., data = train.data)
# make predictions
predictions <- model4 %>% predict(test.data)
# model performance
data.frame(
  RMSE = RMSE(predictions, test.data$W17Height),
  R2 = R2(predictions, test.data$W17Height)
)
model4_summary<-tidy(model4)

# detect multicollinearity
car::vif(model4) 

# detect multicollinearity
VIF <- car::vif(model4) %>%
  as.list() %>% 
  as.data.frame() %>% 
  gather(key = 'variable', value = 'VIF') %>% 
  arrange(desc(VIF))

# remove variables with high VIF (above 5-10)
best.vars<-unique(VIF$variable[which(VIF$VIF<100)])
best.vars

# e.g.
# all vars
SPmod1 <- lmer(W17Height ~ MAT_T + MAP_T + FFP_T + DD5_T + MSP_T + eFFP_T + NFFD_T + MCMT_T + EMNT_T + MWMT_T + AHM_T 
               + (1|Trial) 
               + (1|sample),
               data = pinus)

summary(SPmod1)
coef(SPmod1)
augment(SPmod1)

# plot residuals
plot(SPmod1, which = 1)
# q plot
qqnorm(resid(SPmod1))
qqline(resid(SPmod1))

###################################################
# centre and scale
sp3<-read.csv("./Scots_pine/Scots_pine_H.csv")
sp3$X<-NULL

# remove NAs
sp3<-na.omit(sp3)
head(sp3)
sp3<-sp3[,-c(1:10,12)]
variables<-unique(colnames(sp3))

# using functions from diagnosing_collinearity.R 
for (i in c(1:44)){
  #i<-1
  var1<-variables[i]
  var2<-sp3[, c(var1)] 
  var_cent<-c.(var2)
  var_scale<-z.(var_cent)
  sp3[,c(var1)]<-var_scale
}

sp3<-read.csv("./Scots_pine/Scots_pine_H_cent_scal_allvars.csv")
#sp3 <- read.csv("~/Documents/FR/FR_R/DeltaTrait_bitbucket/Scots_pine/Scots_pine_H_cent_scal_allvars.csv")
summary(sp3)
sp3$X<-NULL
sp3$Latitude<-NULL
sp3$Longitude<-NULL
sp3$Elevation<-NULL

sp <- read.csv("./Scots_pine/Scots_pine_H.csv")
sp <- na.omit(sp)
sp3$Trial <- sp$Trial
sp3$Provenance <- sp$Provenance
sp3$Block <- sp$Block
sp3$Family <- sp$Family
sp3$Seedling <- sp$Seedling
sp3$Tag <- sp$Tag
head(sp3)

corrplot(cor(sp3[,-c(40:46)]), method = "ellipse") 
#pairs(sp3, na.action(na.omit))
collinearity <- as.data.frame(cor(sp3[,-c(40:46)]))

cor <- collinearity %>% 
  mutate(Var1 = factor(row.names(.), levels=row.names(.))) %>% 
  gather(key = Var2, value = value, -Var1, na.rm = TRUE, factor_key = TRUE) 

high.cor <- filter(cor, value <=-0.6 | value >=0.6)
ok.cor <- filter(cor, value >-0.6 & value <=0.6)

#vars<-unique(ok.cor$Var2)
#vars<-as.character(vars)
#sp4 <- sp3[, (names(sp3) %in% vars)]

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

# try performance package
r2(model1)
icc(model1)
check_model(model1)

# detect multicollinearity
car::vif(model1) # there are aliased coeffs 
ld.vars <- attributes(alias(model1)$Complete)$dimnames[[1]]
ld.vars
#  [1] "TD_P"    "FFP_P"   "MCMT_T"  "TD_T"    "MAP_T"   "MSP_T"   "AHM_T"   "SHM_T"   "DD0_T"   "DD5_T"   "DD_18_T" "NFFD_T"  "bFFP_T"  "eFFP_T" 
# [15] "FFP_T"   "PAS_T"   "EMNT_T"  "Eref_T"  "CMD_T"  
sp3<-sp3[,-which(names(sp3) %in% ld.vars)] # remove these

training.samples <- sp3$W17Height %>% createDataPartition(p = 0.8, list = FALSE)
train.data  <- sp3[training.samples, ]
test.data <- sp3[-training.samples, ]

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

# try performance package
r2(model2)
icc(model2)
check_model(model2)

# detect multicollinearity
car::vif(model2)

# detect multicollinearity
VIF <- car::vif(model2) %>%
  as.list() %>% 
  as.data.frame() %>% 
  gather(key = 'variable', value = 'VIF') %>% 
  arrange(desc(VIF))

# remove variables with high VIF (above 5-10)
best.vars<-unique(VIF$variable[which(VIF$VIF<10)])
best.vars

model3 <- lmer(W17Height ~. + (1|Provenance) + (1|Trial/Block/Seedling), data = train.data)

model4 <- lmer(W17Height ~ MWMT_T + (1|Provenance) + (1|Trial/Block/Seedling), data = train.data)
# Random effect variances not available. Returned R2 does not account for random effects.
# Can't compute random effect variances. Some variance components equal zero.
# Solution: Respecify random structure! 
r2(model4)
icc(model4)
check_model(model4)

# detect multicollinearity
car::vif(model4)

# detect multicollinearity
VIF3 <- car::vif(model4) %>%
  as.list() %>% 
  as.data.frame() %>% 
  gather(key = 'variable', value = 'VIF') %>% 
  arrange(desc(VIF))

# remove variables with high VIF (above 5-10)
best.vars<-unique(VIF3$variable[which(VIF3$VIF<10)])
best.vars

### loop to remove each var in turn and see which improves model
#vars<-unique(ok.cor$Var2)
#vars<-as.character(vars)
#vars<-vars[-1]
#model<-c()
#peformance<-c()

#for (i in length(vars)){
#i<-1
#var<-noquote(vars[i])
#model[i] <- lm(W17Height ~. -var, data = train.data)
# model performance
#performance[i] <- model_performance(model[i])
#}

# 09/03/2020
# post-chat with Marta and Richard

sp<-read.csv("./Scots_pine/Scots_pine_H_cent_scal_allvars.csv")
head(sp)
sp$X<-NULL
sp$Latitude<-NULL
sp$Longitude<-NULL
sp$Elevation<-NULL
sp2<- read.csv("./Scots_pine/Scots_pine_H.csv")
sp2<-na.omit(sp2)
sp$Trial <- sp2$Trial
sp$Provenance <- sp2$Provenance
sp$ID <- sp2$ID1
sp$Block <- sp2$Block
sp$Family <- sp2$Family
sp$Seedling <- sp2$Seedling
sp$Tag <- sp2$Tag
sp$row <- sp2$Row
sp$column <- sp2$Column
sp$H2 <- sp2$W17Height
head(sp)
rm(sp2)

# make nesting explicit
sp <- sp %>% mutate(block=factor(Trial:Block),
                    popSite=factor(Trial:Provenance),
                    famSite=factor(Trial:Family))
head(sp)

## overall summary - variation among populations...
ggplot(sp, aes(reorder(Provenance, W17Height), W17Height, group = Provenance))+
  geom_violin()+
  facet_wrap(~Trial)

## spatial effect in trial?
ggplot(sp, aes(row, column, fill = W17Height))+
  geom_tile(colour = "grey55")+theme_void()+
  facet_wrap(~Trial)+scale_fill_viridis_c()

colnames(sp)[1]<-"H"

# boxplots
boxplot(sp$H2 ~ sp$Trial)
par(mar=c(10,4,4,2))
boxplot(sp$H2 ~ sp$ID, las=2)
dotchart(sp$H2, groups=factor(sp$Trial), color = sp$Trial)
dotchart(sp$H2, groups=factor(sp$ID), color = sp$ID)

# relationships between the response variable and the explanatory variables
# pairplots
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

graphics::pairs(sp[,c(1:10)], diag.panel=panel.hist, upper.panel=panel.cor)
graphics::pairs(sp[,c(53,2:10)], diag.panel=panel.hist, upper.panel=panel.cor) # raw data response var

sp$H2log<-log(sp$H2) # log raw height
graphics::pairs(sp[,c(54,2:10)], diag.panel=panel.hist, upper.panel=panel.cor) # log raw height
graphics::pairs(sp[,c(54,11:21)], diag.panel=panel.hist, upper.panel=panel.cor)
graphics::pairs(sp[,c(54,22:31)], diag.panel=panel.hist, upper.panel=panel.cor)
graphics::pairs(sp[,c(54,32:40)], diag.panel=panel.hist, upper.panel=panel.cor)

# start testing combinations of variables

# MAP_T and FFP_P (based on annual precipitation and growing season length being most important components in Richard's thesis analysis)
SPmod1 <- lmer(H ~ MAP_T + FFP_P + (1|Provenance) + (1|Trial/Block), data = sp)
tidy(SPmod1)
r2(SPmod1)
icc(SPmod1)
check_model(SPmod1)

# same but log(H)
SPmod2 <- lmer(log(H) ~ MAP_T + FFP_P + (1|Provenance) + (1|Trial/Block), data = sp)
r2(SPmod2)
icc(SPmod2)
check_model(SPmod2)
summary(SPmod2)
#fixef(SPmod2) # fixed effects
#ranef(SPmod2) # random effects
#confint(SPmod2) # confidence intervals
tidy(SPmod2, conf.int=TRUE) # repeats three previous calls in one call

# switch trial and provenance variables
SPmod3 <- lmer(log(H) ~ MAP_P + FFP_T + (1|Provenance) + (1|Trial/Block), data = sp)
r2(SPmod3)
icc(SPmod3)
check_model(SPmod3)

# all variables from PC1 of Richard's thesis
# growing degree-days, monthly mean temps for Feb and July, annual precipitation, extreme temperature range
SPmod4 <- lmer(log(H) ~ DD5_T + MWMT_P + MCMT_T + MAP_T + TD_P + (1|Provenance) + (1|Trial/Block), data = sp)
model_performance(SPmod4)
check_model(SPmod4)

# comparisons
compare_performance(SPmod1,SPmod2,SPmod3,SPmod4, rank = TRUE)
plot(compare_performance(SPmod1,SPmod2,SPmod3,SPmod4, rank = TRUE))

# variables shown to be correlated with H from pairplots
# o	MWMT_T (0.41)
# o	PAS_T (0.41)
# o	TD_T (0.33)
# o	DD0_T (0.30)
# o	MCMT_T (0.28)
# o	eFFP_T (0.28)
# o	Eref_T (0.24)
# o	EMNT_T (0.22)
# o	DD18_T (0.21)
# o	NFFD_T (0.21) 

# linear mod for all
spLM <- lm(H ~ MWMT_T + PAS_T + TD_T + DD0_T + MCMT_T + eFFP_T + Eref_T + EMNT_T + DD_18_T + NFFD_T, data = sp )
summary(spLM)

# highest two correlated - mean warmest month temp and precipitation as snow
#SPmod5 <- lmer(log(H) ~ MWMT_T + PAS_T + (1|Trial/Provenance/Family/Block), data=sp)
#SPmod5 <- lmer(log(H) ~ MWMT_T + PAS_T + (1|Trial/Block/Provenance/Family), data=sp)
#summary(SPmod5)
#tidy(SPmod5, conf.int=TRUE)

# work out nesting properly
sp <- sp %>% mutate(nest=factor(Trial:Block:Provenance))
length(unique(sp$nest))

# precipitation as snow, and number of frost-free days
SPmod6 <- lmer(log(H) ~ PAS_T + NFFD_T + (1|popSite) + (1|block), data=sp)
#SPmod6a <- lmer(log(H) ~ PAS_T + NFFD_T + (1|Provenance) + (1|Trial/Block), data=sp)
#SPmod6b <- lmer(log(H) ~ PAS_T + NFFD_T + (1|Provenance/Family) + (1|Trial/Block), data=sp)
SPmod6c <- lmer(log(H) ~ PAS_T + NFFD_T + (1|nest), data=sp)

#SPmod7 <- lmer(log(H) ~ PAS_T + NFFD_T + (1|nest), data=sp)

SPmod8 <- lmer(log(H) ~ DD0_T + TD_T + (1|popSite) + (1|block), data=sp)
SPmod9 <- lmer(log(H) ~ PAS_T + DD0_T + TD_T + (1|popSite) + (1|block), data=sp)

compare_performance(SPmod6,SPmod6c,SPmod8, SPmod9, rank = TRUE)
plot(compare_performance(SPmod6,SPmod6c,SPmod8, SPmod9, rank = TRUE))

check_model(SPmod8)

# dredge for all possible combinations
require(MuMIn)
options(na.action = "na.fail") # change the default "na.omit" to prevent models from being fitted to different datasets in case of missing values.

globalmodel <- lmer(log(H2) ~ MWMT_T + PAS_T + TD_T + DD0_T + MCMT_T + eFFP_T + Eref_T + EMNT_T + DD_18_T + NFFD_T +
                    (1|popSite) + (1|block),
                    data = sp)
summary(globalmodel)

combinations <- dredge(globalmodel)
print(combinations)
coefs <- coefTable(combinations)
coefTable(combinations)[1]

combinations<- combinations[order(combinations$AICc),]
models <- get.models(combinations, subset=TRUE)

# "best models"
# growing degree days above 18, number of frost-free days
mod1<-lmer(log(H2) ~ DD_18_T + NFFD_T + (1|popSite) + (1|block),data = sp)
# same and also precipitation as snow
mod2<-lmer(log(H2) ~ DD_18_T + NFFD_T + PAS_T + (1|popSite) + (1|block),data = sp)
# same with temperature difference
mod3<-lmer(log(H2) ~ DD_18_T + NFFD_T + TD_T + (1|popSite) + (1|block),data = sp)
# combination of all of the above
mod4<-lmer(log(H2) ~ DD_18_T + NFFD_T + PAS_T + TD_T + (1|popSite) + (1|block),data = sp)
mod5<-lmer(log(H2) ~ DD_18_T + EMNT_T + (1|popSite) + (1|block),data = sp)
mod6<-lmer(log(H2) ~ DD_18_T + EMNT_T + Eref_T + (1|popSite) + (1|block),data = sp)
mod7<-lmer(log(H2) ~ DD_18_T + EMNT_T + MCMT_T + (1|popSite) + (1|block),data = sp)
mod8<-lmer(log(H2) ~ DD_18_T + EMNT_T + Eref_T + MCMT_T + (1|popSite) + (1|block),data = sp)

# realised EMNT was a typo by me when processing climate EU data
# it's just EMT - extreme minimum temperature

compare_performance(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8, rank = TRUE)
plot(compare_performance(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8, rank = TRUE))

# run on training and test datasets
results<-NULL
for(i in 1:1000) {
  #i<-1
  rows<-sample(x=1:nrow(sp),size=round(nrow(sp)*0.25,1))
  training<-sp[-rows,]
  testing<-sp[rows,]
  mod1<-lmer(log(H2) ~ DD_18_T + NFFD_T + (1|popSite) + (1|block),data = training)
  mod2<-lmer(log(H2) ~ DD_18_T + NFFD_T + PAS_T + (1|popSite) + (1|block),data = training)
  mod3<-lmer(log(H2) ~ DD_18_T + NFFD_T + TD_T + (1|popSite) + (1|block),data = training)
  mod4<-lmer(log(H2) ~ DD_18_T + NFFD_T + PAS_T + TD_T + (1|popSite) + (1|block),data = training)
  mod5<-lmer(log(H2) ~ DD_18_T + EMNT_T + (1|popSite) + (1|block),data = training)
  mod6<-lmer(log(H2) ~ DD_18_T + EMNT_T + Eref_T + (1|popSite) + (1|block),data = training)
  mod7<-lmer(log(H2) ~ DD_18_T + EMNT_T + MCMT_T + (1|popSite) + (1|block),data = training)
  mod8<-lmer(log(H2) ~ DD_18_T + EMNT_T + Eref_T + MCMT_T + (1|popSite) + (1|block),data = training)
  prd1<-predict(mod1,testing)
  errors<-prd1-testing$H
  results<-rbind(results,data.frame(Run=i,MAE=mean(abs(errors)),RMSE=sqrt(sum(errors^2)/nrow(testing))))
}
summary(results)
hist(results$MAE,breaks=50)
hist(results$RMSE,breaks=50)
boxplot(results$MAE)
boxplot(results$RMSE)

# PCA model - try all possible combos
require(MuMIn)
options(na.action = "na.fail") # change the default "na.omit" to prevent models from being fitted to different datasets in case of missing values.

head(sp)
sp$H<-NULL

PCAmodel <- lmer(log(H2) ~ DD_18_T + FFP_T + bFFP_P + FFP_P + TD_P +
                      (1|popSite) + (1|block),
                    data = sp)
summary(PCAmodel)

combinations <- dredge(PCAmodel)
print(combinations)
coefs <- coefTable(combinations)
coefTable(combinations)[1]

combinations<- combinations[order(combinations$AICc),]
models <- get.models(combinations, subset=TRUE)

pcaMOD1<-lmer(log(H2) ~ DD_18_T + FFP_T + (1|popSite) + (1|block),data = sp)
pcaMOD2<-lmer(log(H2) ~ bFFP_P + DD_18_T + FFP_T + (1|popSite) + (1|block),data = sp)

check_model(PCAmodel)
check_model(pcaMOD2)

compare_performance(PCAmodel,pcaMOD1,pcaMOD2, rank = TRUE)
plot(compare_performance(mod1,mod2,mod3,mod4,mod5,mod6,mod7,mod8, rank = TRUE))


# increase understanding of data
# look at height by main explanatory variables
require(Hmisc)
# read in non-standardised data
Pinus<-read.csv("./Scots_pine/Scots_pine_H.csv")
Pinus$X<-NULL
# remove NAs
Pinus<-na.omit(Pinus)
head(Pinus)
summary(Pinus$DD_18_T)
Pinus$DD18_cut<-cut(Pinus$DD_18_T, seq(3600,4000,100))
summary(Pinus$DD18_cut)
dotchart(Pinus$W17Height,groups = factor(Pinus$DD18_cut),color = Pinus$Trial, ylab='DD18 (grouped)', xlab = 'Height', main='Coloured by trial')
dotchart(Pinus$W17Height,groups = factor(Pinus$DD18_cut),color = Pinus$Provenance, ylab='DD18 (grouped)', xlab = 'Height', main='Coloured by provenance')

summary(Pinus$FFP_T)
Pinus$FFP_cut<-cut(Pinus$FFP_T, seq(170,220,10))
summary(Pinus$FFP_cut)
dotchart(Pinus$W17Height,groups = factor(Pinus$FFP_cut),color = Pinus$Trial, ylab='FFP (grouped)', xlab = 'Height', main='Coloured by trial')
dotchart(Pinus$W17Height,groups = factor(Pinus$FFP_cut),color = Pinus$Provenance, ylab='FFP (grouped)', xlab = 'Height', main='Coloured by provenance')

summary(Pinus$bFFP_P)
Pinus$bFFP_cut<-cut(Pinus$bFFP_P, seq(90,150,10))
summary(Pinus$bFFP_cut)
dotchart(Pinus$W17Height,groups = factor(Pinus$bFFP_cut),color = Pinus$Trial, ylab='FFP (grouped)', xlab = 'Height', main='Coloured by trial')
dotchart(Pinus$W17Height,groups = factor(Pinus$bFFP_cut),color = Pinus$Provenance, ylab='FFP (grouped)', xlab = 'Height', main='Coloured by provenance')

summary(Pinus$TD_P)
Pinus$TD_cut<-cut(Pinus$TD_P, seq(10,20,5))
summary(Pinus$TD_cut)
dotchart(Pinus$W17Height,groups = factor(Pinus$TD_cut),color = Pinus$Trial, ylab='FFP (grouped)', xlab = 'Height', main='Coloured by trial')
dotchart(Pinus$W17Height,groups = factor(Pinus$TD_cut),color = Pinus$Provenance, ylab='FFP (grouped)', xlab = 'Height', main='Coloured by provenance')


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
null_model <- lmer(log(H2) ~ (1 | Provenance) , data = sp)

# Build the Trial model with Trial as a fixed and random slope and provenance as the random-effect
trial_model <- lmer(log(H2) ~ Trial + (1 + Trial | Provenance) , data = sp)
block_model <- lmer(log(H2) ~ Block + (1 + Block | popSite), data = sp) # use explicitly nested vars (Block = Trial:Block, popSite=Trial:Provenance)

# Compare null_model and year_model using an anova
anova(null_model, trial_model)
anova(trial_model,block_model)






















