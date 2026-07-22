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
library(effects)
library(car)
library(emmeans)
library(glmmTMB)
library(AICcmodavg)
library(scales)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::summarize)
conflicts_prefer(EnvStats::print)
# set working directory
setwd("~/Documents/EEID_HOFI/git_hofi/code/Katie")

# read data
full <- as_tibble(read.csv("1c_master_20250915.csv"))

full %>% group_by(population) %>% dplyr::summarize(N=n())
full %>% filter(primary_dose!="0") %>% dplyr::summarize(N=n())

seccontrols <-full %>% filter(secondary_dose=="0") 

# subset initial inoculated birds only
full <- full %>% filter(primary_dose!="0") 

# format data
full$secondary_dose <- as.numeric(as.character(full$secondary_dose))
full$primary_dose <- as.numeric(as.character(full$primary_dose))
full$sec_inf <- as.factor(full$sec_inf)
full$sex <- as.factor(full$sex)
full$prim_inf <- as.factor(full$prim_inf)
full$population <- as.factor(full$population)
full$prim_inf <- relevel(full$prim_inf, ref = 2)
levels(full$sex) <- c("F", "M")
levels(full$population) <- c("AZ", "VA")
levels(full$prim_inf) <- c("Infected", "Uninfected")


# specify infection criteria
full$prim_inf <- with(full, ifelse(eyetot7 > 0 | eyetot14 > 0 | eyetot21 > 0 | 
                                     eyetot28 > 0 | eyetot35 > 0 | eyetot41 > 0 |quantity7>15| eyetot41 > 0,  1, 0)) 

full$sec_inf <- with(full, ifelse(eyetot46 > 0 | eyetot56 > 0 | 
                                    quantity46 > 15 | quantity49 > 15 |quantity56>15,  1, 0)) 

full$sec_inf <- with(full, ifelse(band==2426, 0, full$sec_inf)) #fix NA for bird 2426

# calculate initial tolerance
prim_infected <- full %>% filter(prim_inf=="1")
prim_uninfected <- full %>% filter(prim_inf=="0")

presid_mod <- glm(I(prim_toleye/6)~prim_resist10, data=prim_infected, 
                  weights=rep(12,length(prim_infected$prim_resist10)), family="binomial")

prim_infected$presid <- residuals.glm(presid_mod,type="response") #store residuals
prim_infected$presid <- prim_infected$presid*(-1) #flip sign of residuals (now most tolerant residuals are positive)

full <- rbind.fill(prim_infected, prim_uninfected)

# remove second inoculation NAs and controls
full <- full %>% filter(sec_inf!="NA") %>% filter(secondary_dose!="0") 

# calculate second tolerance
sinfected <- full %>% filter(sec_inf==1)#subset primary infected
suninfected <- full %>% filter(sec_inf==0) #subset primary uninfected
suninfected$sresid <- NA

sresid_mod <- glm(I(sec_toleye/6)~sec_resist10, data=sinfected, 
                  weights=rep(12,length(sinfected$sec_resist10)), family="binomial")

sinfected$sresid <- residuals.glm(sresid_mod,type="response") #store residuals
sinfected$sresid <- sinfected$sresid*(-1) #flip sign of residuals (now most tolerant residuals are positive)

full <- rbind.fill(sinfected, suninfected)

# format sec_infected df
sinfected$secondary_dose <- as.numeric(as.character(sinfected$secondary_dose))
sinfected$sex <- as.factor(sinfected$sex)
sinfected$prim_inf <- as.factor(sinfected$prim_inf)
sinfected$population <- as.factor(sinfected$population)

n_table <-full %>% group_by(population, secondary_dose, sex) %>% tally()


# view second infection data
loadnoeye <- full %>% filter(sec_maxeye ==0 & sec_maxload > 0)
min(loadnoeye$sec_maxload) 
  # min load for birds with zero eye score

full %>% filter(prim_maxeye >0 & prim_maxload>0) %>% summarize(minload = min(prim_maxload))
full %>% filter(sec_maxeye >0 & sec_maxload>0) %>% summarize(minload = min(sec_maxload))


full %>% filter(sec_maxeye > 0 & sec_maxload>0) %>% dplyr:: summarize(minload = min(sec_maxload))
  # min load for birds with nonzero eye score and nonzero load
  # note there is one bird with eye score >0 and a load of 0

# List A/Table S2: Probability of second-inoculation infection (all birds) ------------------

# create models
secinf.null <- glm(sec_inf~1, data=full, family="binomial")
secinf.a <- glm(sec_inf~secondary_dose, data=full, family="binomial") 

secinf.b <- glm(sec_inf~secondary_dose+population, data=full, family="binomial")
secinf.c <- glm(sec_inf~secondary_dose*population, data=full, family="binomial")

secinf.d <- glm(sec_inf~secondary_dose+sex, data=full, family="binomial") 
secinf.e <- glm(sec_inf~secondary_dose*sex, data=full, family="binomial") 

secinf.f <- glm(sec_inf~secondary_dose+prim_inf, data=full, family="binomial")
secinf.g <- glm(sec_inf~secondary_dose*prim_inf, data=full, family="binomial")

secinf.h <- glm(sec_inf~secondary_dose+prim_resist10, data=full, family="binomial")
secinf.i <- glm(sec_inf~secondary_dose*prim_resist10, data=full, family="binomial")

secinf.j <- glm(sec_inf~secondary_dose+prim_toleye, data=full, family="binomial")
secinf.k <- glm(sec_inf~secondary_dose*prim_toleye, data=full, family="binomial")

secinf.l <- glm(sec_inf~secondary_dose+population+sex, data=full, family="binomial") 
secinf.m <- glm(sec_inf~secondary_dose+population*sex, data=full, family="binomial")

secinf.n <- glm(sec_inf~secondary_dose+prim_inf*population, data=full, family="binomial")
secinf.o <- glm(sec_inf~secondary_dose+prim_inf+population, data=full, family="binomial")

secinf.p <- glm(sec_inf~secondary_dose+prim_inf*sex, data=full, family="binomial")
secinf.q <- glm(sec_inf~secondary_dose+prim_inf+sex, data=full, family="binomial")

secinf.r <- glm(sec_inf~secondary_dose+prim_toleye*sex, data=full, family="binomial")
secinf.s <- glm(sec_inf~secondary_dose+prim_toleye+sex, data=full, family="binomial")

secinf.t <- glm(sec_inf~secondary_dose+prim_toleye*population, data=full, family="binomial")
secinf.u <- glm(sec_inf~secondary_dose+prim_toleye+population, data=full, family="binomial")

secinf.v <- glm(sec_inf~secondary_dose+prim_resist10*population, data=full, family="binomial")
secinf.w <- glm(sec_inf~secondary_dose+prim_resist10+population, data=full, family="binomial")

secinf.x <- glm(sec_inf~secondary_dose+prim_resist10*sex, data=full, family="binomial")
secinf.y <- glm(sec_inf~secondary_dose+prim_resist10+sex, data=full, family="binomial")


# model list
secinf_list <- list(secinf.null=secinf.null, secinf.a=secinf.a, secinf.b=secinf.b, secinf.c=secinf.c, 
                    secinf.d=secinf.d, secinf.e=secinf.e, secinf.f=secinf.f, secinf.g=secinf.g,
                    secinf.h=secinf.h, secinf.i=secinf.i, secinf.j=secinf.j, secinf.k=secinf.k,
                    secinf.l=secinf.l, secinf.m=secinf.m, secinf.n=secinf.n, secinf.o=secinf.o,
                    secinf.p=secinf.p, secinf.q=secinf.q, secinf.r=secinf.r, secinf.s=secinf.s,
                    secinf.t=secinf.t, secinf.u=secinf.u, secinf.v=secinf.v, secinf.w=secinf.w,
                    secinf.x=secinf.x, secinf.y=secinf.y)

aictab(cand.set=secinf_list) 

# check out models <2 dAICc
summary(secinf.d) # top model
anova(secinf.d, test="Chisq")
Anova(secinf.d)
res <- simulateResiduals(secinf.d)
plot(res) 

summary(secinf.e) # more complex than top model (interaction term)
anova(secinf.e, test="Chisq")
Anova(secinf.e)
res <- simulateResiduals(secinf.e)
plot(res) 

summary(secinf.l) 
anova(secinf.l, test="Chisq") # p<0.10
Anova(secinf.l)
res <- simulateResiduals(secinf.l)
plot(res) 

summary(secinf.a) # nested in top model
anova(secinf.a, test="Chisq")
Anova(secinf.a)
res <- simulateResiduals(secinf.a)
plot(res) 

summary(secinf.y)
anova(secinf.y, test="Chisq") # p<0.10
Anova(secinf.y)
res <- simulateResiduals(secinf.y)
plot(res) 

summary(secinf.null)


# check relationship btw initial and second susceptibility by chi square 
resp_table = table(full$prim_inf, full$sec_inf) 
EnvStats::print(resp_table)
print(chisq.test(resp_table))



# List B/Table S3: Probability of second-inoculation infection (initial infected birds)  --------

# create models
infected <- full %>% filter(prim_inf=="1") 

secinf2.null <- glm(sec_inf~1, data=infected, family="binomial")
secinf2.a <- glm(sec_inf~secondary_dose, data=infected, family="binomial") 

secinf2.b <- glm(sec_inf~secondary_dose+population, data=infected, family="binomial")
secinf2.c <- glm(sec_inf~secondary_dose*population, data=infected, family="binomial")

secinf2.d <- glm(sec_inf~secondary_dose+sex, data=infected, family="binomial") 
secinf2.e <- glm(sec_inf~secondary_dose*sex, data=infected, family="binomial") 

secinf2.f <- glm(sec_inf~secondary_dose+prim_resist10, data=infected, family="binomial")
secinf2.g <- glm(sec_inf~secondary_dose*prim_resist10, data=infected, family="binomial")

secinf2.h <- glm(sec_inf~secondary_dose+prim_toleye, data=infected, family="binomial")
secinf2.i <- glm(sec_inf~secondary_dose*prim_toleye, data=infected, family="binomial")

secinf2.j <- glm(sec_inf~secondary_dose+population+sex, data=infected, family="binomial")  
secinf2.k <- glm(sec_inf~secondary_dose+population*sex, data=infected, family="binomial")

secinf2.l <- glm(sec_inf~secondary_dose+prim_toleye*sex, data=infected, family="binomial")
secinf2.m <- glm(sec_inf~secondary_dose+prim_toleye+sex, data=infected, family="binomial")
secinf2.n <- glm(sec_inf~secondary_dose+prim_toleye*population, data=infected, family="binomial")
secinf2.o <- glm(sec_inf~secondary_dose+prim_toleye+population, data=infected, family="binomial")

secinf2.p <- glm(sec_inf~secondary_dose+prim_resist10*population, data=infected, family="binomial")
secinf2.q <- glm(sec_inf~secondary_dose+prim_resist10+population, data=infected, family="binomial")
secinf2.r <- glm(sec_inf~secondary_dose+prim_resist10*sex, data=infected, family="binomial")
secinf2.s <- glm(sec_inf~secondary_dose+prim_resist10+sex, data=infected, family="binomial")

secinf2.t <- glm(sec_inf~secondary_dose+presid, data=infected, family="binomial")
secinf2.u <- glm(sec_inf~secondary_dose+presid+population, data=infected, family="binomial")
secinf2.v <- glm(sec_inf~secondary_dose+presid+sex, data=infected, family="binomial") 
secinf2.w <- glm(sec_inf~secondary_dose+presid+sex+population, data=infected, family="binomial") 

# model list
secinf_list2 <- list(secinf2.null=secinf2.null, secinf2.a=secinf2.a, secinf2.b=secinf2.b, secinf2.c=secinf2.c, 
                    secinf2.d=secinf2.d, secinf2.e=secinf2.e, secinf2.f=secinf2.f, secinf2.g=secinf2.g,
                    secinf2.h=secinf2.h, secinf2.i=secinf2.i, secinf2.j=secinf2.j, secinf2.k=secinf2.k,
                    secinf2.l=secinf2.l, secinf2.m=secinf2.m, secinf2.n=secinf2.n, secinf2.o=secinf2.o,
                    secinf2.p=secinf2.p, secinf2.q=secinf2.q, secinf2.r=secinf2.r, secinf2.s=secinf2.s,
                    secinf2.t=secinf2.t, secinf2.u=secinf2.u, secinf2.v=secinf2.v, secinf2.w=secinf2.w)

# model comparison
aictab(cand.set=secinf_list2)

# check out models <2 dAICc
summary(secinf2.w) # top model
anova(secinf2.w, test="Chisq")
Anova(secinf2.w)
secinf2.w_res <- simulateResiduals(secinf2.w)
plot(secinf2.w_res) 

summary(secinf2.j) # nested in top model
anova(secinf2.j, test="Chisq")
Anova(secinf2.j)
secinf2.j_res <- simulateResiduals(secinf2.j)
plot(secinf2.j_res) 

summary(secinf2.v) # nested in top model
anova(secinf2.v, test="Chisq")
Anova(secinf2.v)
secinf2.v_res <- simulateResiduals(secinf2.v)
plot(secinf2.v_res) 

summary(secinf2.null)



# List C/Table S4: Second infection resistance (all second infected birds) --------

# write models
secresist.null <- lm(sec_resist10~1, data=sinfected)

secresist.a <- lm(sec_resist10~secondary_dose, data=sinfected) 
secresist.b <- lm(sec_resist10~population, data=sinfected)
secresist.c <- lm(sec_resist10~sex, data=sinfected)
secresist.d <- lm(sec_resist10~prim_inf, data=sinfected)
secresist.e <- lm(sec_resist10~prim_resist10, data=sinfected)
secresist.f <- lm(sec_resist10~prim_toleye, data=sinfected)

secresist.g <- lm(sec_resist10~secondary_dose+population, data=sinfected) 
secresist.h <- lm(sec_resist10~secondary_dose+sex, data=sinfected) 
secresist.i <- lm(sec_resist10~prim_toleye+population, data=sinfected) 
secresist.j <- lm(sec_resist10~secondary_dose+prim_inf, data=sinfected)
secresist.k<- lm(sec_resist10~secondary_dose+prim_resist10, data=sinfected)
secresist.l<- lm(sec_resist10~secondary_dose+prim_toleye, data=sinfected)
secresist.m <- lm(sec_resist10~prim_inf+population, data=sinfected) 
secresist.n <- lm(sec_resist10~prim_inf+sex, data=sinfected) 
secresist.o <- lm(sec_resist10~prim_resist10+population, data=sinfected) 
secresist.p <- lm(sec_resist10~prim_resist10+sex, data=sinfected) 
secresist.q<- lm(sec_resist10~prim_toleye+sex, data=sinfected) 

secresist.r <- lm(sec_resist10~sex+population+secondary_dose+prim_resist10, data=sinfected)
secresist.s<- lm(sec_resist10~sex+population+secondary_dose+prim_toleye, data=sinfected)
secresist.t <- lm(sec_resist10~sex+population+secondary_dose+prim_inf, data=sinfected)


# model list
secresist_list <- list(secresist.null=secresist.null, secresist.a=secresist.a, secresist.b=secresist.b, secresist.c=secresist.c, 
                       secresist.d=secresist.d, secresist.e=secresist.e, secresist.f=secresist.f,
                       secresist.g=secresist.g, secresist.h=secresist.h, secresist.i=secresist.i, secresist.j=secresist.j,
                       secresist.k=secresist.k,
                       secresist.l=secresist.l, secresist.m=secresist.m, secresist.n=secresist.n, secresist.o=secresist.o,
                       secresist.p=secresist.p, secresist.q=secresist.q, secresist.r=secresist.r, secresist.s=secresist.s,
                       secresist.t=secresist.t)

# model comparison
aictab(cand.set=secresist_list) 

# check out model <2 dAICc
summary(secresist.j) # top and only model <2 dAICc
#plot(secresist.j)

summary(secresist.null)



# List D/Table S5: Second infection resistance (second infected AZ birds)  ------------

sinfected_az <- sinfected %>% filter(population=="AZ") 

# write models
secresist.null2 <- lm(sec_resist10~1, data=sinfected_az)

secresist.a2 <- lm(sec_resist10~secondary_dose, data=sinfected_az) 
secresist.b2 <- lm(sec_resist10~sex, data=sinfected_az)
secresist.c2 <- lm(sec_resist10~prim_inf, data=sinfected_az)
secresist.d2 <- lm(sec_resist10~prim_resist10, data=sinfected_az)
secresist.e2 <- lm(sec_resist10~prim_toleye, data=sinfected_az)


# model list
secresist_list2 <- list(secresist.null2=secresist.null2, secresist.a2=secresist.a2, 
                       secresist.b2=secresist.b2, secresist.c2=secresist.c2, 
                       secresist.d2=secresist.d2, secresist.e2=secresist.e2)

# model comparison
aictab(cand.set=secresist_list2)

# check out models <2 dAICc
summary(secresist.c2)
#plot(secresist.c2)

summary(secresist.d2)
#plot(secresist.d2)

summary(secresist.e2)
#plot(secresist.e2)

summary(secresist.null2)

# List E/Table S6: Second infection tolerance (all second infected birds) --------

# create tolerance residuals 
sresid_mod <- glm(I(sec_toleye/6)~sec_resist10, data=sinfected, 
                 weights=rep(12,length(sinfected$sec_resist10)), family="binomial")
summary(sresid_mod)

sinfected$sresid <- residuals.glm(sresid_mod,type="response") 
sinfected$sresid <- sinfected$sresid*(-1) 
  #flip sign of residuals (now most tolerant residuals are positive)

# write models
sectol.null <- lm(sresid~1, data=sinfected) 

sectol.a <- lm(sresid~secondary_dose, data=sinfected) 
sectol.b <- lm(sresid~population, data=sinfected)
sectol.c <- lm(sresid~sex, data=sinfected)
sectol.d <- lm(sresid~prim_inf, data=sinfected)
sectol.e <- lm(sresid~prim_resist10, data=sinfected)
sectol.f <- lm(sresid~prim_toleye, data=sinfected)

sectol.g <- lm(sresid~secondary_dose+population, data=sinfected) 
sectol.h <- lm(sresid~secondary_dose+sex, data=sinfected) 
sectol.i <- lm(sresid~prim_toleye+population, data=sinfected) 
sectol.j <- lm(sresid~secondary_dose+prim_inf, data=sinfected)
sectol.k<- lm(sresid~secondary_dose+prim_resist10, data=sinfected)
sectol.l<- lm(sresid~secondary_dose+prim_toleye, data=sinfected)

sectol.m <- lm(sresid~prim_inf+population, data=sinfected) 
sectol.n <- lm(sresid~prim_inf+sex, data=sinfected) 
sectol.o <- lm(sresid~prim_resist10+population, data=sinfected) 
sectol.p <- lm(sresid~prim_resist10+sex, data=sinfected) 
sectol.q<- lm(sresid~prim_toleye+sex, data=sinfected) 
sectol.r <- lm(sresid~sex+population+secondary_dose+prim_resist10, data=sinfected)
sectol.s<- lm(sresid~sex+population+secondary_dose+prim_toleye, data=sinfected)
sectol.t <- lm(sresid~sex+population+secondary_dose+prim_inf, data=sinfected)

# model list
sectol_list <- list(sectol.null=sectol.null, sectol.a=sectol.a, sectol.b=sectol.b, sectol.c=sectol.c, 
                       sectol.d=sectol.d, sectol.e=sectol.e, sectol.f=sectol.f,
                       sectol.g=sectol.g, sectol.h=sectol.h, sectol.i=sectol.i, sectol.j=sectol.j,
                       sectol.k=sectol.k,
                       sectol.l=sectol.l, sectol.m=sectol.m, sectol.n=sectol.n, sectol.o=sectol.o,
                       sectol.p=sectol.p, sectol.q=sectol.q, sectol.r=sectol.r, sectol.s=sectol.s)

# model comparison
aictab(cand.set=sectol_list)
  # null is top ranked by weight

summary(sectol.null) # "top" model

# summaries for Table S6
summary(sectol.a)
summary(sectol.b)
summary(sectol.d)




# List F/Table S7: Second infection tolerance (second infected AZ birds) ------------------

sinfected_az <- sinfected %>% filter(population=="AZ") 

# write models
sectol2.null <- lm(sresid~1, data=sinfected_az)

sectol2.a <- lm(sresid~secondary_dose, data=sinfected_az) 
sectol2.b <- lm(sresid~sex, data=sinfected_az)
sectol2.c <- lm(sresid~prim_inf, data=sinfected_az)
sectol2.d <- lm(sresid~prim_resist10, data=sinfected_az)
sectol2.e <- lm(sresid~prim_toleye, data=sinfected_az)
sectol2.f <- lm(sresid~secondary_dose+sex, data=sinfected_az)
sectol2.g <- lm(sresid~secondary_dose+prim_inf, data=sinfected_az)
sectol2.h<- lm(sresid~secondary_dose+prim_resist10, data=sinfected_az)
sectol2.i<- lm(sresid~secondary_dose+prim_toleye, data=sinfected_az)
sectol2.j <- lm(sresid~prim_inf+sex, data=sinfected_az) 
sectol2.k <- lm(sresid~prim_resist10+sex, data=sinfected_az) 
sectol2.l<- lm(sresid~prim_toleye+sex, data=sinfected_az) 

# model list
sectol_list2 <- list(sectol2.null=sectol2.null, sectol2.a=sectol2.a,  sectol2.b=sectol2.b, 
                     sectol2.c=sectol2.c, sectol2.d=sectol2.d, sectol2.e=sectol2.e, sectol2.f=sectol2.f,
                     sectol2.g=sectol2.g,sectol2.h=sectol2.h, sectol2.i=sectol2.i, sectol2.j=sectol2.j,
                     sectol2.k=sectol2.k, sectol2.l=sectol2.l)

# model comparison
aictab(cand.set=sectol_list2)

# check out model <2 dAICc
summary(sectol2.g) # top model
#plot(sectol2.g)

summary(sectol2.c) # nested in top model

summary(sectol2.null)





# Figure 3 ----------------------------------------------------------------

# Figure 3A: second infection probability by host sex and second dose (all birds)
colors <- c("#648FFF", "#FFB000")
pinfected <- full %>% filter(prim_inf == "1")

par(mfrow=c(2,1), mai=c(1,1,0.5,0.2))
# mai =c(bottom, left, top, right)
visreg(secinf.d, xvar = "secondary_dose", by="sex", scale="response", partial=F, 
       overlay=TRUE, cex.axis=0.75,   ylim=c(-0.1, 1.1), xtrans=log10,
       xlab="Second MG dose (CCU/mL)", ylab="Second infection probability", rug=F,
       line=list(col=c("#648FFF", "#FFB000")),
       fill=list(col=c("#648FFF80", "#FFB00080"))) #, xaxt="n")
title("A", adj=0)
axis(1, at=c(1.477121, 2, 2.477121, 3.845098),
     labels=c("30", "100", "300", "7,000"))
points(jitter(sec_inf, 0.4) ~ jitter(log10(secondary_dose), 0.8), full, pch = 1, 
       bg = 'grey', cex = 0.8, col= colors[factor(sex)])
# doses are 30 (1.477121), 100(2), 300 (2.477121), and 7000 CCU (3.845098)






# Figure 3B: second infection probability by initial tolerance (initial infected birds only)
visreg(secinf2.w, "presid", scale="response", partial=F,
       overlay=FALSE, cex.axis=0.75, ylim = c(-0.1,1.1), rug=F,
       xlab="Initial tolerance", ylab="Second infection probability",
       line=list(col=c("black")))
points(jitter(sec_inf, 0.2) ~ presid, pinfected, pch = 1, bg = 'grey', cex = 0.8)
title("B", adj=0)


# Figure 4 ----------------------------------------------------------------

group.colors <- c( "#785EF0","#FE6100")
pred_data <- expand.grid(secondary_dose = seq(from = min(sinfected$secondary_dose, na.rm = T), 
                                              to = max(sinfected$secondary_dose, na.rm = T), by = 30),
                         prim_inf = factor(c('0', '1'), levels = c('0', '1')),
                         sec_resist10 = seq(from = min(sinfected$sec_resist10, na.rm = T), 
                                            to = max(sinfected$sec_resist10, na.rm = T), by = 0.05))

cis<-stats::predict(secresist.j, newdata = pred_data, type='response', se.fit=TRUE) 
pred_data$pred <- cis$fit
pred_data$se <- cis$se.fit
pred_data$upr.ci <- pred_data$pred + 1.96*pred_data$se
pred_data$lwr.ci <- pred_data$pred - 1.96*pred_data$se

ggplot(sinfected, aes(x = secondary_dose, y = sec_resist10, fill=prim_inf)) +
  geom_jitter(width=0.05, height=0, size=2, aes(shape=factor(population), color=prim_inf)) + 
  scale_x_continuous(trans="log10", breaks=c(30,100, 300,7000))+
  scale_y_continuous(breaks=c(0,1,3,5), labels=c("0", "10", "1000", "100000")) +
  geom_line(data = pred_data, aes(x = secondary_dose, y = pred, col = prim_inf))+
  labs(y= "Second MG load + 1", x = "Second MG dose (CCU/mL)", 
       color="Initial \nsusceptibility") +
  scale_color_manual(values=group.colors, labels=c("Uninfected", "Infected")) +
  theme(legend.position = "right", 
        panel.border=element_rect(colour="black", fill=NA, linewidth=1.5),
        panel.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.text=element_text(size=12))+
  geom_ribbon(data=pred_data, aes(ymin=lwr.ci, ymax=upr.ci), 
              alpha=0.3, linetype=0)+
  scale_fill_manual(values=group.colors, guide="none")+
  scale_shape(solid=FALSE)+
  labs(shape="Host \npopulation", colour="Initial \nsusceptibility")


# Figure 5 ----------------------------------------------------------------

labels<-c("Unsusceptible", "Susceptible")

ggplot(sinfected, aes(x = prim_inf, y = sresid)) +
  labs(y= "Second inoculation tolerance", x = "Initial susceptibility", 
       shape="Host \npopulation")+
  theme(legend.position = "right",
        panel.border=element_rect(colour="black", fill=NA, linewidth=1.5),
        panel.background=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        legend.title = element_text(size=14),
        legend.text=element_text(size=12))+
  geom_jitter(width=0.09, height=0, size=3, aes(shape=factor(population)))+
  scale_shape(solid=FALSE) + scale_x_discrete(label=labels) 




# Compare initial and second inoculation severity -------------------------

# how many finches infected during both inoculations?
full %>% filter(prim_inf=="1" & sec_inf=="1") %>% summarize(n=n())

full %>% group_by(prim_inf, sec_inf) %>%summarize(n=n(), primload.m=mean(prim_resist10), primload.se =(sd(prim_resist10))/(sqrt(n)),
                        secload.m=mean(sec_resist10), secload.se=(sd(sec_resist10))/(sqrt(n)),
            primeye.m=mean(prim_maxeye), primeye.se=(sd(prim_maxeye))/(sqrt(n)),
            seceye.m=mean(sec_maxeye), seceye.se=(sd(sec_maxeye))/(sqrt(n)))

full %>% filter(prim_inf=="1" & sec_inf=="0") %>% group_by(secondary_dose) %>% summarize(n=n())
full %>% filter(prim_inf=="0" & sec_inf=="1") %>% group_by(secondary_dose) %>% summarize(n=n())
full %>% filter(prim_inf=="1" & sec_inf=="1") %>% group_by(secondary_dose) %>% summarize(n=n())


infecteds <- full %>% filter(prim_inf=="1" & sec_inf=="1")
prim_resp <- infecteds[c(1,2,4,34,35,36, 43)]
prim_resp <- prim_resp %>% filter(prim_infected==1)
sec_resp <- infecteds[c(1,2,4,40,41,42, 44)]
prim_resp$inoc <- "primary"
sec_resp$inoc <- "secondary"

prim_resp <-dplyr::rename(prim_resp, maxload=prim_resist10, maxeye=prim_toleye,infstatus=prim_inf, tol=presid)
sec_resp <- dplyr::rename(sec_resp, maxload=sec_resist10, maxeye=sec_toleye, infstatus=sec_inf, tol=sresid)

sevdata <- rbind(prim_resp, sec_resp)

hist(sevdata$maxload)
hist(sevdata$maxeye)
hist(sevdata$tol)

sev_load <- lmer(maxload~inoc + (1|band), data=sevdata)
sev_load


#Fig 1 style primary to secondary inoc resp
ggplot(full, aes(x = prim_resist10, y = sec_resist10)) +
  labs(y= "Second inoculation resistance", x = "Initial resistance")+
  geom_jitter(width=0.04, height=0.04)

ggplot(infecteds, aes(x=presid, y=sresid))+
  geom_point()

ggplot(full, aes(x=prim_maxeye, y=sec_maxeye))+
  geom_jitter(height=0.04, width=0.04)

ggplot(full, aes(x=presid, y=sec_resist10))+
  geom_jitter(height=0.04, width=0.04)
