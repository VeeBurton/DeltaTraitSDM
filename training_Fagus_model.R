
# date: 17/02/20
# author: VB
# description: using Fagus climate and trait data, learn how to check for colinearity and choose the best variables for the model
# will need to do this for poplar and scots pine
# see what i can do with Fagus and compare to Marta's out the box model in Height_2070_RCP8.5.R

# libraries
library(tidyverse)

# read in data
Fagus <- read.table("./Fagus/Fagus_H_Final.txt",header=T)
str(Fagus)


