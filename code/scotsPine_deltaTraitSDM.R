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

# based on all analysis so far, these are the most important variables according to:

# pairplots:
    # MWMT_T (0.41)
    # PAS_T (0.41)
    # TD_T (0.33)
    # DD0_T (0.30)
    # MCMT_T (0.28)
    # eFFP_T (0.28)
    # Eref_T (0.24)
    # EMT_T (0.22)
    # DD18_T (0.21)
    # NFFD_T (0.21) 
# PCA: 
    # DD_18_T
    # FFP_T
    # bFFP_P
    # FFP_P
    # TD_P
# exploration/difference between trial & provenance:
    # PAS
    # TD
    # MWMT
    # DD0
    # MCMT

# read in standardised data (with raw height)
sp.raw <- read.csv("./Scots_pine/Scots_pine_H.csv") # raw data
sp.raw$X <- NULL
str(sp.raw)
sp <- read.csv("./Scots_pine/Scots_pine_H_cent_scal_allvars.csv") # 
sp$X <- NULL
str(sp)
colnames(sp)[1]<-"H"
sp<-na.omit(sp)
sp$Hraw <- sp.raw$H
sp$Trial <- sp.raw$Trial
sp$Provenance <- sp.raw$Provenance
sp$Block <- sp.raw$Block

