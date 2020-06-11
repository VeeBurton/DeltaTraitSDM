library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(lme4)
library(corrplot)

sp <- read.csv("~/R/DeltaTraitSDM/Scots_pine/PS_Common_Garden_Heights.csv")

colnames(sp)[3]<-"Trial"
colnames(sp)[8]<-"Provenance"
colnames(sp)[11]<-"height"
sp$ï..Seqno<-NULL
summary(sp$height)
sp$height <- as.numeric(as.character(sp$height)) 
summary(sp$height)
# i had erronous height values before due to incorrect conversion from factor to numeric
sp$Block[which(sp$Block=="c")] <- "C"
sp$Block<-factor(sp$Block)

summary(sp)

# fully understand structure
sp.structure <- sp %>% 
  group_by(Trial,Provenance) %>% 
  summarise(Provenances=length(unique(Provenance)),
            Families=length(unique(as.character(Family))),
            Individuals=length(unique(Tag)))

sp.families <- sp %>% 
  group_by(Provenance) %>% 
  summarise(Families=list(unique(Family)))
# definitely unique Families per provenance

# make sure random effect stucture is correctly defined
# data is crossed - we observed every Provenance and Family at every Trial and Block
# there are 3 sites – called Borders, Glensaugh, Inverewe here – and each is a fully randomised block design: 
# there are 4 blocks at Borders, Glensaugh and 3 at Inverewe
# in each site there are 21 populations, 8 families and either 3 or 4 individuals per family (1 per block)
# so there are 21*8*4=672 at Borders, Glensaugh; 21*8*3= 504 at Inverewe; total = 1848

# eplicitly nest blocks within trials
sp <- within(sp, site <- factor(Trial:Block))
length(unique(sp$site)) # 11, 4 @ Borders, 4 @ Glensaugh, 3 @ Inverewe
# could represent nesting either with (1|Trial/Block) or (1|site)
xtabs(~ Trial + Block, sp) # incorrect
xtabs(~ Trial + site, sp) # correct
# explcitly nest families within provenances
sp <- within(sp, provFam <- factor(Provenance:Family))
length(unique(sp$provFam)) # mostly 8 per provenance (21*8=168), some have 9 families
xtabs(~ Provenance + Family, sp)
xtabs(~ Provenance + provFam, sp)

# every provenance and every family is at every trial and block (fully crossed)
xtabs(~ Trial + provFam, sp)
xtabs(~ site + provFam, sp)
site.mod <- lmer(height ~ (1|site)+(1|provFam), sp) # crossed design with explicitly nested variables
trial.mod <- lmer(height ~ (1|Trial)+(1|provFam), sp) # crossed design ignoring block, ProvFam explicitly nested
summary(site.mod)
summary(trial.mod)
# site (Trial:Block) is used, there one observation of ProvFam per site
# if Trial is used, there are 3-4 observations of ProvFam per site

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
sp2$Trial <- factor(sp2$Trial)
sp2$Provenance <- factor(sp2$Provenance)

summary(sp2$height)
hist(sp2$height)
# normal distribution
summary(sp2)

write.csv(sp2,"~/R/DeltaTraitSDM/Scots_pine/sp_data_height_climateTP.csv")
sp2 <- read.csv("~/R/DeltaTraitSDM/Scots_pine/sp_data_height_climateTP.csv")
sp2<-na.omit(sp2)

# look at height distribution per provenance for each trial
ggplot(sp2, aes(Provenance, height, color = Provenance))+
  geom_violin()+
  facet_grid(Trial~Block)+
  theme(axis.text.x = element_text(angle = 90))+
  theme(legend.position = "none")

# family means
sp.summary <- sp2 %>% na.omit %>%  
  group_by(Trial,provFam) %>% 
  summarise(mn_height=mean(height)) 

sp.summary2 <- sp.summary %>%
  group_by(Trial) %>% 
  arrange(desc(mn_height), .by_group = TRUE)

write.csv(sp.summary2, "~/R/DeltaTraitSDM/Scots_pine/sp_family_summaries.csv")

# plot family means per trial
ggplot(sp.summary)+
  geom_col(aes(provFam,mn_height,fill=provFam))+
  facet_wrap(~Trial)+
  theme(legend.position = "none")

# two way anova for significant differences between family means
bxp <- ggboxplot(
  sp2, x = "Trial", y = "height",
  color = "provFam")
bxp
# check normality by examining residuals
# Build the linear model
lmod  <- lm(height ~ Trial*provFam, data = sp2)
summary(lmod)
# Create a QQ plot of residuals
ggqqplot(residuals(lmod)) # along line, normal
ggqqplot(sp, "height", ggtheme = theme_bw()) +
  facet_grid(provFam ~ Trial)
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(lmod)) # if the p-value is not significant, we can assume normality
# homogneity of variance assumption
# this can be checked using the Levene’s test, if not significant assume homogeneity
sp2 %>% levene_test(height ~ Trial*provFam) # not significant
res.aov <- sp2 %>% anova_test(height ~ Trial*provFam)
res.aov
# statistically significant interaction between Trial and provFam (F=1.2), p<0.05*

# simple mixed model
mxmod <- lmer(height ~ (1|Trial)+(1|provFam), sp)
summary(mxmod)
# according to (Trial) + (Provenance/Family)
# Trial explains 51% of the variance
# Provenance:Family explains 6.8% of the variance
# there is still 41% residual error

# plot observed height vs. fitted height
sp2 %>% cbind(fitted = fitted(mxmod))%>%
  group_by(Trial, provFam, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%
  filter(diff>=500 | diff<=-500) %>% 
  ggplot(aes(Trial, diff, group = provFam, colour = Seed.Zone, label = provFam))+
  geom_line()+geom_text()+
  ylab("Residuals")

sp2%>%cbind(fitted = fitted(mxmod))%>%
  group_by(Trial, Provenance, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%#filter(Trial != "BORDERS")%>%
  tidyr::pivot_wider(id_cols = c("Provenance", "Seed.Zone"),
                     names_from = Trial,
                     values_from = diff)%>%
  mutate(sum_diff = abs(GLENSAUGH) + abs(INVEREWE))%>%
  filter(sum_diff > 300)%>%
  reshape2::melt(id.vars = c("Provenance", "Seed.Zone", "sum_diff"))%>%
  ggplot(aes(variable, value, group = Provenance, colour = Seed.Zone, label = Provenance))+
  geom_line()+geom_text()+geom_hline(yintercept = 0, lty = "dashed")+
  ylab("Residuals")+xlab("Trial")

# residuals = fitted - observed
# positive residuals mean height prediction was too low
# negative residuals mean height prediction was too high

sp2%>%cbind(fitted = fitted(mxmod))%>%
  filter(Trial!="GLENSAUGH") %>% 
  ggplot(aes(fitted,height, color=Trial))+geom_point()+geom_smooth()

residFit <- sp2 %>% na.omit %>% cbind(fitted = fitted(mxmod))%>%
  group_by(Trial, Provenance, Seed.Zone)%>%
  summarise(mean_fitted = mean(na.omit(fitted)),
            mean_obs = mean(na.omit(height)))%>%
  mutate(diff = mean_fitted-mean_obs)%>%
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

# explore simple relationships 
# e.g.
ggplot(sp2, 
       aes(x = cut(DD0_P, breaks = 5), y = height, colour=Provenance)) + 
  geom_boxplot()#geom_point()

# sometimes there really is no meaningful relationship between the two variables. 
# other times, a careful transformation of one or both of the variables can reveal a clear relationship.
sp2$Latitude.y<-NULL
sp2$Longitude.y<-NULL
sp2$Elevation.y<-NULL
sp_long <- sp2 %>% 
  pivot_longer(cols = MAT_P:CMD_T,
               names_to = "variables",
               values_to = "values")
# reaction norms
ggplot(sp_long, aes(x = values, y = height, color=Trial)) +
  geom_point() +
  geom_smooth(method="lm")+
  facet_wrap(~variables, scales="free_x")+
  scale_x_log10() + 
  scale_y_log10("Height")+
  theme_bw()

# mean differences between provenances and trials
sp3 <- sp2 %>% mutate(MATdiff = MAT_T-MAT_P,
                      MWMTdiff = MWMT_T-MWMT_P,
                      MCMTdiff = MCMT_T-MCMT_P,
                      TDdiff = TD_T-TD_P,
                      MAPdiff = MAP_T-MAP_P,
                      MSPdiff = MSP_T-MSP_P,
                      AHMdiff = AHM_T-AHM_P,
                      SHMdiff = SHM_T-SHM_P,
                      DD0diff = DD0_T-DD0_P,
                      DD_18diff = DD_18_T-DD_18_P,
                      NFFDdiff = NFFD_T-NFFD_P,
                      bFFPdiff = bFFP_T-bFFP_P,
                      eFFPdiff = eFFP_T-eFFP_P,
                      FFPdiff = FFP_T-FFP_P,
                      PASdiff = PAS_T-PAS_P,
                      EMTdiff = EMT_T-EMT_P,
                      Erefdiff = Eref_T-Eref_P,
                      CMDdiff = CMD_T-CMD_P)
mean_diffs <- sp3 %>% 
  group_by(Trial,Provenance) %>% 
  summarise(meanH = mean(height, na.rm=TRUE),
            MAP = mean(MAPdiff),
            MAT = mean(MATdiff),
            TD = mean(TDdiff),
            DD0 = mean(DD0diff),
            DD_18 = mean(DD_18diff),
            NFFD = mean(NFFDdiff),
            bFFP = mean(bFFPdiff),
            FFP = mean(FFPdiff),
            eFFP = mean(eFFPdiff),
            PAS = mean(PASdiff),
            Eref = mean(Erefdiff),
            CMD = mean(CMDdiff))

mean_diffs <- mean_diffs %>% 
  pivot_longer(cols = MAP:CMD,
               names_to = "variables",
               values_to = "meandiff")
mean_diffs$variables <- factor(mean_diffs$variables)
mean_diffs$variables <- ordered(mean_diffs$variables, levels = c("MAP","MAT","TD", "DD0" ,"DD_18",
                                                                 "NFFD","bFFP","FFP","eFFP","PAS",
                                                                 "Eref","CMD"))

ggplot(mean_diffs, aes(meandiff, meanH, colour = Trial))+
  geom_point()+geom_smooth()+facet_wrap(~variables, scales = "free")+
  theme_bw()+theme(legend.position = "bottom")+
  xlab("Climate distance (Trial climate - Provenance climate)")+
  ylab("Mean height (mm)")+
  geom_smooth(data = mean_diffs, aes(x = meandiff, y = meanH, group = 1), lty = "dashed", method = "lm", colour = "black", se = F)

# height ~ all climate variables
# centre and scale climate variables 
sp2$Latitude.y<-NULL
sp2$Longitude.y<-NULL
sp2$Elevation.y<-NULL
sp2[,c(17:56)] <- scale(sp2[,c(17:56)], scale = TRUE, center = TRUE)
summary(sp2)
sp2$DD18_T<-NULL #(mostly NA values)

# corrplots
corrplot(cor(sp2[,c(10,17:36)]), method = "ellipse") # provenance climate
corrplot(cor(sp2[,c(10,37:55)]), method = "ellipse") # trial climate

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

graphics::pairs(sp2[,c(10,17:26)], diag.panel=panel.hist, upper.panel=panel.cor)
graphics::pairs(sp2[,c(10,27:36)], diag.panel=panel.hist, upper.panel=panel.cor)
graphics::pairs(sp2[,c(10,37:46)], diag.panel=panel.hist, upper.panel=panel.cor)
graphics::pairs(sp2[,c(10,47:55)], diag.panel=panel.hist, upper.panel=panel.cor)


mean_diffs <- sp3 %>% 
  group_by(Trial,Provenance) %>% 
  summarise(meanH = mean(height, na.rm=TRUE),
            MAP = mean(MAPdiff),
            MAT = mean(MSPdiff),
            TD = mean(TDdiff),
            DD0 = mean(DD0diff),
            DD18 = mean(DD_18diff),
            NFFD = mean(NFFDdiff),
            bFFP = mean(bFFPdiff),
            FFP = mean(FFPdiff),
            eFFP = mean(eFFPdiff),
            PAS = mean(PASdiff),
            Eref = mean(Erefdiff),
            CMD = mean(CMDdiff))

mean_diffs <- mean_diffs %>% 
  pivot_longer(cols = MAP:CMD,
               names_to = "variables",
               values_to = "meandiff")

ggplot(mean_diffs, aes(meandiff, meanH, colour = Trial))+
  geom_point()+geom_smooth()+facet_wrap(~variables, scales = "free")+
  theme_bw()+theme(legend.position = "bottom")+
  geom_smooth(data = mean_diffs, aes(x = meandiff, y = meanH, group = 1), lty = "dashed", method = "lm", colour = "black", se = F)

# DTSDM

mod1 <- lmer(height ~ MAP_T + NFFD_T + MAP_T*MAP_P + NFFD_T*NFFD_P + (1|Trial)+(1|provFam), sp2)
summary(mod1)
