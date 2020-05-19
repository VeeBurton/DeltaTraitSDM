
install.packages("parameters")
library("parameters")
library(tidyverse)

sp <- read.csv(paste0(wd,"data-raw/Scots_pine_H_cent_scal_allvars.csv"))
str(sp)
names(sp)

# test
lm(Hraw ~ TD_P*DD5_T*NFFD_T, data = sp) %>% 
  select_parameters() %>% 
  model_parameters()
