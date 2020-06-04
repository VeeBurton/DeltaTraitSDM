library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(lme4)

sp <- read.csv("~/R/DeltaTraitSDM/Scots_pine/PS_Common_Garden_Heights.csv")

colnames(sp)[3]<-"Trial"
colnames(sp)[8]<-"Provenance"
colnames(sp)[11]<-"height"
sp$ï..Seqno<-NULL
summary(sp$height)
sp$height <- as.numeric(as.character(sp$height)) 
# i had erronous height values before due to incorrect conversion from factor to numeric
sp$Block[which(sp$Block=="c")] <- "C"
sp$Block<-factor(sp$Block)

summary(sp)

env <- read.csv("~/R/DeltaTraitSDM/Scots_pine/scots_pine_all_locations_elev_Normal_1961_1990Y.csv")
env_P <- env %>% filter(ID2!="Trial")
colnames(env_P)[2]<-"Provenance"
env_T <- env %>% filter(ID2 == "Trial")
colnames(env_T)[1]<-"Trial"
env_T$Trial<-as.character(env_T$Trial)
env_T$Trial[which(env_T$Trial=="Yair")]<-"BORDERS"
env_T$Trial[which(env_T$Trial=="Inverewe")]<-"INVEREWE"
env_T$Trial[which(env_T$Trial=="Glensaugh")]<-"GLENSAUGH"
env_T$ID2<-NULL

colnames(env_P)[6:25] <- c('MAT_P','MWMT_P','MCMT_P','TD_P','MAP_P','MSP_P','AHM_P','SHM_P','DD0_P','DD5_P','DD_18_P','DD18_P','NFFD_P','bFFP_P','eFFP_P','FFP_P','PAS_P','EMT_P','Eref_P','CMD_P')
colnames(env_T)[5:24] <- c('MAT_T','MWMT_T','MCMT_T','TD_T','MAP_T','MSP_T','AHM_T','SHM_T','DD0_T','DD5_T','DD_18_T','DD18_T','NFFD_T','bFFP_T','eFFP_T','FFP_T','PAS_T','EMT_T','Eref_T','CMD_T')

sp2 <- left_join(sp, env_P)
sp2 <- right_join(sp2, env_T, by="Trial")

summary(sp2)

summary(sp2$height)
hist(sp2$height)
# normal distribution

# look at height distribution per provenance for each trial
ggplot(sp2, aes(Provenance, height, color = Provenance))+
  geom_violin()+
  facet_grid(Trial~Block)+
  theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position = "none")

# family means
sp.summary <- sp2 %>% na.omit %>%  
  group_by(Trial,Provenance,Family) %>% 
  summarise(mn_height=mean(height)) 

sp.summary <- sp.summary %>%
  group_by(Trial) %>% 
  arrange(desc(mn_height), .by_group = TRUE)

write.csv(sp.summary, "~/R/DeltaTraitSDM/Scots_pine/sp_family_summaries.csv")

# plot family means per trial
ggplot(sp.summary)+
  geom_col(aes(Provenance,mn_height,fill=Provenance))+
  facet_wrap(~Trial)+
  theme(legend.position = "none")

# three-way anova for significant differences between family means
# check normality by examining residuals
# Build the linear model
lmod  <- lm(height ~ Trial*Provenance*Family, data = sp2)
summary(lmod)
# Create a QQ plot of residuals
ggqqplot(residuals(lmod)) # along line, normal
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(lmod)) # if the p-value is not significant, we can assume normality
# homogneity of variance assumption
# this can be checked using the Levene’s test, if not significant assume homogeneity
sp2 %>% levene_test(height ~ Trial*Provenance*Family)
res.aov <- sp2 %>% anova_test(height ~ Trial*Provenance*Family)
res.aov

# simple mixed model
mxmod <- lmer(height ~ (1|Trial/Provenance/Family), sp2)
summary(mxmod)

total.var <- sum(94517,6306,9371,70808)

Trial<- 94517/total.var*100
Trial.prov<- 6306/total.var*100
Trial.prov.fam <- 9371/total.var*100
Residual<-70808/total.var*100

# trial explains most of the variation (52%)
# family explains more variation (5%) than Provenance (3.5%)
# still a lot of residual (39%)

# plot observed height vs. fitted height
sp2 %>% na.omit %>% cbind(fitted = fitted(mxmod))%>%
  group_by(Trial, Provenance, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%
  ggplot(aes(Trial, diff, group = Provenance, colour = Seed.Zone, label = Provenance))+
  geom_line()+geom_text()

## there is clearly some interaction going on - just look at Glensaugh and Inverewe and get 
## rid of the less interactive/ill-fitting provenances (abritrary cutoff of 25 mm across both sites)

sp2 %>% na.omit %>% cbind(fitted = fitted(mxmod))%>%
  group_by(Trial, Provenance, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%filter(Trial != "BORDERS")%>%
  tidyr::pivot_wider(id_cols = c("Provenance", "Seed.Zone"),
                     names_from = Trial,
                     values_from = diff)%>%
  mutate(sum_diff = abs(GLENSAUGH) + abs(INVEREWE))%>%
  filter(sum_diff > 25)%>%
  reshape2::melt(id.vars = c("Provenance", "Seed.Zone", "sum_diff"))%>%
  ggplot(aes(variable, value, group = Provenance, colour = Seed.Zone, label = Provenance))+
  geom_line()+geom_text()+geom_hline(yintercept = 0, lty = "dashed")

residFit <- sp2 %>% na.omit %>% cbind(fitted = fitted(mxmod))%>%
  group_by(Trial, Provenance, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%filter(Trial != "BORDERS")%>%
  tidyr::pivot_wider(id_cols = c("Provenance", "Seed.Zone"),
                     names_from = Trial,
                     values_from = diff)%>%
  reshape2::melt(id.vars = c("Provenance", "Seed.Zone")) %>% 
  left_join(sp2)

ggplot(residFit, aes(Longitude.x, value, colour = variable, label = Provenance))+
  geom_smooth()+
  geom_smooth(method = "lm", lty = "dashed", se = F)+
  geom_label()+
  theme_bw(base_size = 8)+labs(x = "Longitude", y = "Fitted-Observed")

# height ~ all climate variables

basic.lm <- lm(height ~ Trial:Provenance:Family, data = sp2)
summary(basic.lm)
# plot the data
ggplot(sp3, aes(x =Trial:Provenance, y = H)) +
  geom_jitter() +
  geom_smooth(method = "lm")
# plot residuals
plot(basic.lm, which = 1) # red line should be flat like the dashed grey line
# look at qqplot
plot(basic.lm, which = 2) # points should fall on diagonal line

# check for data independence - essentially checking if need to account for random effects
boxplot(Hraw ~ Trial:Provenance, data = sp3)
# colour plot
ggplot(sp3, aes(x = Trial:Provenance, y = Hraw, colour = Trial:Provenance)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(legend.position = "none")
# split plot
ggplot(aes(Trial:Provenance, Hraw), data = sp3) + 
  geom_point() + 
  facet_wrap(~ Trial:Provenance) + # create a facet for each trial
  xlab("Trial:Provenance") + 
  ylab("height")

# include provenance as a fixed effect
prov.lm <- lm(Hraw ~ Trial:Provenance, data = sp3)
tidy(prov.lm)
# as a random effect
mixed.lmer <- lmer(Hraw ~ DD_18_T + (1|Trial/Provenance), data = sp3)
summary(mixed.lmer)
tidy(mixed.lmer)
plot(mixed.lmer) 
qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer)) # points should fall on line



