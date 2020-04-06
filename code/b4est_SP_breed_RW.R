##Scots pine B4EST

library(dplyr)
library(ggplot2)

#read in data and call interaction/nested variables....

psy = read.csv("Y:/CONIFERS/SCOTS PINE/psy_Scav.csv", na.strings = "*")%>%
  mutate(block = factor(PlantingSite:Block),
         PopSite = factor(PlantingSite:Population),
         FamSite = factor(PlantingSite:Family))

# read in the environmental covariates and join them to the phenotype data

env = read.csv("Y:/CONIFERS/SCOTS PINE/sp_env.csv")

psy = left_join(psy, env)

## overall summary - variation among populations...

ggplot(psy, aes(reorder(Population, W17Height), W17Height, group = Population))+
  geom_violin()+
  facet_wrap(~PlantingSite)

## spatial effect in trial?

ggplot(psy, aes(Row, Column, fill = W17Height))+
  geom_tile(colour = "grey55")+theme_void()+
  facet_wrap(~PlantingSite)+scale_fill_viridis_c()

## melt data to plot response against each covariate...

psySum = psy %>%
  tidyr::pivot_longer(cols = x:VAPPressure)%>%
  group_by(Population, PlantingSite, name, value)%>%
  summarise(mn_ht = mean(na.omit(W17Height)))
  
ggplot(psySum, aes(value, mn_ht, colour = PlantingSite))+
  geom_point()+geom_smooth()+facet_wrap(~name, scales = "free")+
  theme_bw()+theme(legend.position = "bottom")+
  geom_smooth(data = psySum, aes(x = value, y = mn_ht, group = 1), lty = "dashed", method = "lm", colour = "black", se = F)


###

### breedR analysis

library(breedR)

#build pedigree... (for use later)

ped = data.frame(
  dam = as.numeric(psy$Family),
  sire = rep(0, nrow(psy)))
  
ped = ped%>%
  cbind(self = as.numeric(seq(from = max(ped$dam)+1, to = nrow(ped)+max(ped$dam))))%>%
  select(self, sire, dam)

psyDat = psy%>%
  cbind(ped)

##Start by fitting a simple model with just provenance and planting site...

resNest <- remlf90(fixed  = W17Height ~ 1,
                   random = ~ Population + PlantingSite + Block,
                   data = psyDat)

##Shows variance components - still a lot of residual but PlantingSite
##clearly contains most of the variation

summary(resNest)

## Here we join the fitted values and work out residuals associated with each pop and plot them...

psyDat%>%cbind(fitted = fitted(resNest))%>%
  group_by(PlantingSite, Population, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(W17Height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%
  ggplot(aes(PlantingSite, diff, group = Population, colour = Seed.Zone, label = Population))+
  geom_line()+geom_text()

## There is clearly some interaction going on - let's just look at Glensaugh and Inverewe and get 
## rid of the lesss interactive/ill-fitting provenances (abritrary cutoff of 95 mm across both sites)

psyDat%>%cbind(fitted = fitted(resNest))%>%
  group_by(PlantingSite, Population, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(W17Height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%filter(PlantingSite != "BORDERS")%>%
  tidyr::pivot_wider(id_cols = c("Population", "Seed.Zone"),
                     names_from = PlantingSite,
                     values_from = diff)%>%
  mutate(sum_diff = abs(GLENSAUGH) + abs(INVEREWE))%>%
  filter(sum_diff > 90)%>%
  reshape2::melt(id.vars = c("Population", "Seed.Zone", "sum_diff"))%>%
  ggplot(aes(variable, value, group = Population, colour = Seed.Zone, label = Population))+
  geom_line()+geom_text()+geom_hline(yintercept = 0, lty = "dashed")

## So Allt Cul, Glen Affric and Abernethy so better than expected by the model at Glensaugh,
## Beinn Eighe and Shieldaig do better than expected at Inverewe
## Amat does poorer at both.. must be because it did well at Yair...

residFit = psyDat%>%cbind(fitted = fitted(resNest))%>%
  group_by(PlantingSite, Population, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(W17Height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%filter(PlantingSite != "BORDERS")%>%
  tidyr::pivot_wider(id_cols = c("Population", "Seed.Zone"),
                     names_from = PlantingSite,
                     values_from = diff)%>%
  reshape2::melt(id.vars = c("Population", "Seed.Zone"))%>%
  left_join(env)

ggplot(residFit, aes(x, value, colour = variable, label = Population))+
  geom_smooth()+
  geom_smooth(method = "lm", lty = "dashed", se = F)+
  geom_label()+
  theme_bw(base_size = 8)+labs(x = "Longitude OSGB36", y = "Fitted-Observed")
                           
## tendency is for the ones that grow more than predicted by the model at Inverewe are from the west and v.versa
## Glen Affric and Rothie looks a bit funny - plot again without it - just for laughs

residFit%>%filter(Population != "GA" & Population != "RM")%>%
  ggplot(aes(x, value, colour = variable, label = Population))+
  geom_point()+
  geom_smooth()+geom_smooth(method = "lm", lty = "dashed", se = F)

####
###
####

psyDat%>%group_by(Population, PlantingSite, x)%>%
  summarise(mn = mean(na.omit(W17Height)))%>%
  ggplot(aes(x, mn, colour = PlantingSite))+
  geom_point()+geom_smooth(method = "lm")



### site by site breedR with spatial... First Yair

##ignore all for now

#yair = psyDat%>%filter(PlantingSite == "BORDERS")

#gen.psy <- list(model    = 'add_animal', 
#                pedigree = yair[,c("self", "sire", "dam")],
#                id       = 'self')

#resSpGen <- remlf90(fixed  = W17Height ~ 1,
#                      random = ~ Population,
#                      genetic = list(model = 'competition',
#                                pedigree = yair[,c("self", "sire", "dam")],
#                                id = 'self',
#                                coord = yair[, c('Row', 'Column')],
#                                competition_decay = 1,
#                                pec = list(present = TRUE)), 
#                      spatial = list(model = 'AR', 
#                                   coord = yair[, c('Row','Column')]), 
#                      data = yair, method = "em")

#summary(resSpGen)

#yair%>%cbind(fitted = fitted(resSpGen))%>%
#  ggplot(aes(reorder(Population, fitted), W17Height))+
#  geom_boxplot()
