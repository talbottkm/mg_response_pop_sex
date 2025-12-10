# Prepare -----------------------------------------------------------------

# clean off workspace
rm(list=ls()) #clear workspace
cat("\014") # clear console

# library
library(dplyr)
library(readxl)
library(tibble)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(conflicted)
library(tidyverse)
library(boot)
library(DHARMa)
library(visreg)
library(AICcmodavg)
library(ggpubr)
library(glm.predict)
library(EnvStats)
library(emmeans)
library(plyr)
library(dplyr)
library(lme4)
library(MASS)
library(ggplot2)
#library(pscl)
library(effects)
library(car)
library(emmeans)
library(glmmTMB)
library(bbmle)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)

# set working directory
setwd("~/Documents/EEID_HOFI/git_hofi/code/Katie")

# read data
master <- as_tibble(read.csv("1c_master_20250915.csv")) 

# format data
master$primary_dose <- as.numeric(as.character(master$primary_dose))
master$sex <- as.factor(master$sex)
master$prim_inf <- as.factor(master$prim_inf)
master$population <- as.factor(master$population)
master$prim_resist10 <- as.numeric(master$prim_resist10)
master$prim_resist10 <- round(master$prim_resist10, digits=0)

master <- master %>% filter(primary_dose=="750")

master$prim_inf <- with(master, ifelse(eyetot7 > 0 | eyetot14 > 0 | eyetot21 > 0 | 
                                     eyetot28 > 0 | eyetot35 > 0 | eyetot41 > 0 |quantity7>15| eyetot41 > 0,  1, 0)) 

pinfected <- master %>% filter(prim_inf=="1") 
pinfected$group <- paste(pinfected$sex, pinfected$population)
pinfected$group <- as.factor(pinfected$group)

resids <- function(a) {
  b <- simulateResiduals(fittedModel = a, n=1000, plot=TRUE)
  plot(b, asFactor=F)}

levels(master$sex) <- c("F", "M")
levels(master$population) <- c("AZ", "VA")
levels(pinfected$group) <-c("F AZ", "F VA", "M AZ", "M VA")

view <- master %>% dplyr::group_by(primary_dose, secondary_dose, population, sex) %>% dplyr::summarize(n=n())

#calculate initial tolerance 
presid_mod <- glm(I(prim_toleye/6)~prim_resist10, data=pinfected, 
                  weights=rep(12,length(pinfected$prim_resist10)), family="binomial")

pinfected$presid2 <- residuals.glm(presid_mod,type="response") #store residuals
pinfected$presid <- pinfected$presid2*(-1) #flip sign of residuals (now most tolerant residuals are positive)

# Model comparison: Initial inoculation susceptibility ------------------

# create models
priminf.null <- glm(prim_inf~1, data=master, family="binomial")
priminf.a <- glm(prim_inf~population, data=master, family="binomial")
priminf.b <- glm(prim_inf~sex, data=master, family="binomial")
priminf.c <- glm(prim_inf~population+sex, data=master, family="binomial")
priminf.d <- glm(prim_inf~population*sex, data=master, family="binomial")

# model list
priminf_list <- list(priminf.null=priminf.null, priminf.a=priminf.a, priminf.b=priminf.b, priminf.c=priminf.c, 
                     priminf.d=priminf.d)

# model comparison
aictab(cand.set=priminf_list)
  # null is top model

summary(priminf.null)
summary(priminf.b)
summary(priminf.a)

priminf.b_res <- simulateResiduals(priminf.b)
plot(priminf.b_res) 

priminf.a_res <- simulateResiduals(priminf.a)
plot(priminf.a_res) 

# Model comparison: Initial resistance (infected birds only) ----------------------------------

pinfected$prim_resist10 <- round(pinfected$prim_resist10, 0)
pinfected$quantity7 <- round(pinfected$quantity7, 0)

hist(pinfected$prim_resist10)
hist(pinfected$quantity7)

# create models
primresist.null <- glmmTMB(quantity7~1, data=pinfected, family=nbinom1)
primresist.a <- glmmTMB(quantity7~population, data=pinfected, family=nbinom1)
primresist.b <-glmmTMB(quantity7~sex, data=pinfected, family=nbinom1)
primresist.c <- glmmTMB(quantity7~sex+population, data=pinfected, family=nbinom1)
primresist.d <- glmmTMB(quantity7~sex*population, data=pinfected, family=nbinom1)

# model list
primresist_list <- list(primresist.null=primresist.null, primresist.a=primresist.a, 
                        primresist.b=primresist.b, primresist.c=primresist.c, 
                        primresist.d=primresist.d)

# model comparison
aictab(cand.set=primresist_list)
  # null is in top model set; model a is < 2 dAIC

summary(primresist.a)
summary(primresist.null)
summary(primresist.d)
summary(primresist.c)


res = simulateResiduals(primresist.a)
plot(res, rank = T)

res = simulateResiduals(primresist.d)
plot(res, rank = T)

res = simulateResiduals(primresist.c)
plot(res, rank = T)

# try kruskall-wallis to compare resistance by group
kw <-kruskal.test(quantity7~group, data=pinfected)
kw # X2(3)=11.04, p=0.01
FSA::dunnTest(quantity7~group, data=pinfected, method="bh")
  # AZ females vary from VA females (Z=3.25, BH adjusted p=0.007)

# Figure S1
options(scipen = 999)
ggplot(pinfected, aes(x = group, y = prim_resist)) +
  geom_jitter(width=0.1, height=0) +theme(legend.position = "right")+
  labs(x="Group (Sex and population)", y="Initial MG load + 1")+
  scale_y_continuous(labels=scales::comma)

# Model comparison: Initial pathology (ordinal regression) ----------------------------------------------

# includes only infected birds

# model list
glmo0<-polr(as.factor(prim_toleye)~1, data= pinfected, Hess=TRUE,model=TRUE)    
glmo1<-polr(as.factor(prim_toleye)~prim_resist10, data= pinfected, Hess=TRUE,model=TRUE)  
glmo2<-polr(as.factor(prim_toleye)~prim_resist10+population, data= pinfected, Hess=TRUE,model=TRUE)    
glmo3<- polr(as.factor(prim_toleye)~prim_resist10*population, data=pinfected, Hess=TRUE, model=TRUE)
glmo4<- polr(as.factor(prim_toleye)~prim_resist10+sex, data=pinfected, Hess=TRUE, model=TRUE)
glmo5<- polr(as.factor(prim_toleye)~prim_resist10*sex, data=pinfected, Hess=TRUE, model=TRUE)
glmo6<- polr(as.factor(prim_toleye)~prim_resist10+sex+population, data=pinfected, Hess=TRUE, model=TRUE)
glmo7<- polr(as.factor(prim_toleye)~prim_resist10+sex*population, data=pinfected, Hess=TRUE, model=TRUE)
glmo8<- polr(as.factor(prim_toleye)~prim_resist10*sex+population, data=pinfected, Hess=TRUE, model=TRUE)

# model list
primpath_list <- list(glmo0=glmo0, glmo1=glmo1, glmo2=glmo2, glmo3=glmo3,
                      glmo4=glmo4, glmo5=glmo5, glmo6=glmo6, glmo7=glmo7,
                      glmo8=glmo8)

# model comparison
aictab(cand.set=primpath_list)

# summary
summary(glmo2)
summary(glmo3)
summary(glmo0)

# Figure S2
group.colors <- c("#648FFF", "#FFB000")
  # blue female, orange male

ggplot(pinfected, aes(x = population, y = prim_toleye, col=sex)) +
  geom_jitter(width=0.1, height=0.04) +theme(legend.position = "right")+
  labs(x="Population", y="Initial eyescore", title="Initial eye score by population and sex",
       color="Sex")+
  scale_color_manual(values=group.colors, labels=c("F", "M")) +
  scale_x_discrete(labels = c('AZ \n(MG-naive)', 'VA \n(MG-endemic)'))

# Model comparison: Initial tolerance -----------------------------------

# create models
primtol.null <- lm(presid~1, data=pinfected)

primtol.a <- lm(presid~population, data=pinfected)
primtol.b <- lm(presid~sex, data=pinfected)
primtol.c <- lm(presid~population*sex, data=pinfected)
primtol.d <- lm(presid~population+sex, data=pinfected)

# model list
primtol_list <- list(primtol.null=primtol.null, primtol.a=primtol.a, primtol.b=primtol.b, primtol.c=primtol.c, 
                     primtol.d=primtol.d)

# model comparison
aictab(cand.set=primtol_list)

summary(primtol.a)
#plot(primtol.a)

summary(primtol.d)
summary(primtol.null)

# Figure S3
group.colors <- c("#648FFF", "#FFB000")
  # blue female, orange male

ggplot(pinfected, aes(x = population, y = presid, col=sex)) +
  geom_jitter(width=0.1, height=0) +theme(legend.position = "right")+
  labs(x="Population", y="Initial tolerance", title = "Initial tolerance by population and sex",
       color="Sex") +
  scale_color_manual(values=group.colors,labels=c("F", "M")) +
  scale_x_discrete(labels = c('AZ \n(MG-naive)', 'VA \n(MG-endemic)'))









