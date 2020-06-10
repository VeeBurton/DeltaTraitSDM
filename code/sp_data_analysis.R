
sp <- read.csv("~/R/DeltaTraitSDM/Scots_pine/sp_data_height_climateTP.csv")
head(sp)
sp<-na.omit(sp)

# order trials by longitude
sp$Trial <- ordered(sp$Trial, levels = c("INVEREWE", "BORDERS", "GLENSAUGH"))
# provenance boxplot by longitude
ggplot(sp)+
  geom_boxplot(aes(Longitude,height,color=ID1))+
  facet_wrap(~Trial)

# Perry et al. 2016
# To test for significant differences in height among families, populations and blocks
# nested analysis of variance (ANOVA) tests were performed with population as a fixed effect, 
# and families nested within population and block as random effects 

mod1 <- lmer(height ~ Provenance + (1|Trial/Block/Provenance), sp)
summary(mod1)

# George et al. 2020
# fitted a general linear model between height and distance between trial site climate and provenance climate
# transfer distance = prov climate - trial climate (or vice versa?)
