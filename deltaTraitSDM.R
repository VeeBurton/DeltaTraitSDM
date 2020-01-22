# libraries
###################
library(lme4)
library(MuMIn)
library(raster)
library(rgdal)
library(rworldmap)
library(sjPlot) 
library(curry)
library(tidyverse)
library(broom)
###################

# climatic and trait database for Fagus
Fagus <- read.csv("~/R/DeltaTraitSDM-Code/Fagus_DeltaTraitSDM_H.csv",header=T)
head(Fagus)

# Model that splits the climate of the site/climate of the provenance/interaction between both as fixed effects 
# and keeps the provenance and trial/block/tree as random effects

M1 <- lmer(log(H) ~ St_Age # height as a function of stand age
           + St_Bio5_P # and max temp of warmest month
           + St_Pet.max_T  # potential evapotranspiration - max?
           + I(St_Bio5_P^2) # I() prevents conversion to factors
           + I(St_Pet.max_T^2)  
           + St_Age*St_Bio5_P 
           + St_Age*St_Pet.max_T 
           + St_Bio5_P*St_Pet.max_T 
           + (1|Trial/Block/Tree_ID) # and trial location - RANDOM EFFECT
           + (1|ID_ProvCode), # and provenance - RANDOM EFFECT
           Fagus) # data

summary(M1)

# with no log
M1_nolog <- lmer(H ~ St_Age 
                 + St_Bio5_P 
                 + St_Pet.max_T 
                 + I(St_Bio5_P^2) 
                 + I(St_Pet.max_T^2)  
                 + St_Age*St_Bio5_P 
                 + St_Age*St_Pet.max_T 
                 +  St_Bio5_P*St_Pet.max_T 
                 + (1|Trial/Block/Tree_ID) 
                 + (1|ID_ProvCode), 
                 Fagus)

summary(M1_nolog)

#M20 = St_Pet.max_P ~ St_Bio13_T
M20 <- lmer(log(H) ~ St_Age 
            + St_Pet.max_P 
            + St_Bio13_T 
            + I(St_Pet.max_P^2) 
            + I(St_Bio13_T^2)  
            + St_Age*St_Pet.max_P 
            + St_Age*St_Bio13_T 
            +  St_Pet.max_P*St_Bio13_T 
            + (1|Trial/Block/Tree_ID) 
            + (1|ID_ProvCode),
            Fagus)

save(M20, file = "M20.rda")


#Select the model with the lowest AIC
AIC(M1, M20, M1_nolog)

### Cross-validation with independent data - to assess capacity of extrapolation of the model

calibrate.data<-sample(1:nrow(Fagus), 0.66*nrow(Fagus))

Fagus1=Fagus[calibrate.data,]

Fagus2=Fagus[-calibrate.data,]

M1 <- lmer(log(H) ~ St_Age 
           + St_Pet.max_P 
           + St_Bio13_T 
           + I(St_Pet.max_P^2) 
           + I(St_Bio13_T^2)  
           + St_Age*St_Pet.max_P 
           + St_Age*St_Bio13_T 
           + St_Pet.max_P*St_Bio13_T 
           + (1|Trial/Block/Tree_ID) # random effects distinguished by vertical bars |
           + (1|ID_ProvCode), 
           Fagus1)

M2 <- predict (M1, 
              newdata=Fagus2,re.form=NA, type="response")

cor.test(Fagus2$H, M2, method = ("pearson"))


