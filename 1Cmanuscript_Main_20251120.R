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

# set working directory
setwd("~/Documents/EEID_HOFI/git_hofi/code/Katie")

# read data
full <- as_tibble(read.csv("1c_master_20250915.csv"))

full %>% group_by(population) %>% dplyr::summarize(N=n())
full %>% filter(primary_dose!="0") %>% dplyr::summarize(N=n())

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



# Second inoculation susceptibility ------------------

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

secinf.z <- glm(sec_inf~secondary_dose+prim_inf:presid, data=full, family="binomial")
secinf.aa <- glm(sec_inf~secondary_dose+prim_inf:presid+sex, data=full, family="binomial")
secinf.bb <- glm(sec_inf~secondary_dose+prim_inf:presid+population, data=full, family="binomial")

secinf.full <- glm(sec_inf~secondary_dose+sex+population+prim_inf+prim_resist10+prim_toleye+prim_inf:presid, data=full, family="binomial") 


# model list
secinf_list <- list(secinf.null=secinf.null, secinf.a=secinf.a, secinf.b=secinf.b, secinf.c=secinf.c, 
                    secinf.d=secinf.d, secinf.e=secinf.e, secinf.f=secinf.f, secinf.g=secinf.g,
                    secinf.h=secinf.h, secinf.i=secinf.i, secinf.j=secinf.j, secinf.k=secinf.k,
                    secinf.l=secinf.l, secinf.m=secinf.m, secinf.n=secinf.n, secinf.o=secinf.o,
                    secinf.p=secinf.p, secinf.q=secinf.q, secinf.r=secinf.r, secinf.s=secinf.s,
                    secinf.t=secinf.t, secinf.u=secinf.u, secinf.v=secinf.v, secinf.w=secinf.w,
                    secinf.x=secinf.x, secinf.y=secinf.y, secinf.full=secinf.full,
                    secinf.z=secinf.z, secinf.aa=secinf.aa, secinf.bb=secinf.bb)

aictab(cand.set=secinf_list)

# check out top models
summary(secinf.aa)
anova(secinf.aa, test="Chisq")
Anova(secinf.aa)
secinf.aa_res <- simulateResiduals(secinf.aa)
plot(secinf.aa_res) 

summary(secinf.z)
anova(secinf.z, test="Chisq")
secinf.z_res <- simulateResiduals(secinf.z)
plot(secinf.z_res)

summary(secinf.bb)
anova(secinf.bb, test="Chisq")
secinf.bb_res <- simulateResiduals(secinf.bb)
plot(secinf.bb_res)

summary(secinf.null)

# check relationship btw initial and second susceptibility by chi square 
resp_table = table(full$prim_inf, full$sec_inf) 
print(resp_table)
print(chisq.test(resp_table))


# Figure 3A
colors <- c("#648FFF", "#FFB000")
pinfected <- full %>% filter(prim_inf == "1")

par(mfrow=c(1,1))
  # mai =c(bottom, left, top, right); mai=c(1,1,0.5,0.2)
visreg(secinf.aa, "secondary_dose", by="sex", scale="response", partial=F, 
       overlay=TRUE, cex.axis=0.75,  xtrans=log10, ylim=c(-0.1, 1.1),
       xlab="Second MG dose (CCU/mL)", ylab="Second infection probability", rug=F,
       line=list(col=c("#648FFF", "#FFB000")),
       fill=list(col=c("#648FFF80", "#FFB00080")), xaxt="n")
axis(1, at=c(1.477121, 2, 2.477121, 3.845098),
     labels=c("30", "100", "300", "7,000"))
points(jitter(sec_inf, 0.4) ~ jitter(log10(secondary_dose), 0.8), full, pch = 1, 
       bg = 'grey', cex = 0.8, col=colors)
  # doses are 30 (1.477121), 100(2), 300 (2.477121), and 7000 CCU (3.845098)


# Figure 3B
par(mfrow=c(1,1))
visreg(secinf.aa, "presid", scale="response", partial=F,
       overlay=FALSE, cex.axis=0.75, ylim = c(-0.1,1.1), rug=F, cond=list(prim_inf=1),
       xlab="Initial tolerance", ylab="Second infection probability")
points(jitter(sec_inf, 0.2) ~ presid, pinfected, pch = 16, bg = 'grey', cex = 0.8)
  #above uses 'pinfected' to overlay points, and 'full' to create confidence interval 





# Second inoculation susceptibility (initial infected)  --------

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

# check out top model
summary(secinf2.w)
anova(secinf2.w, test="Chisq")
Anova(secinf2.w)
secinf2.w_res <- simulateResiduals(secinf2.w)
plot(secinf2.w_res) 

summary(secinf2.j)
anova(secinf2.j, test="Chisq")
secinf2.j_res <- simulateResiduals(secinf2.j)
plot(secinf2.j_res) 

summary(secinf2.v)
anova(secinf2.v, test="Chisq")
secinf2.v_res <- simulateResiduals(secinf2.v)
plot(secinf2.v_res) 

summary(secinf2.null)

# Second infection resistance (second infected birds) --------

# additional formatting
full$secondary_dose <- as.numeric(as.character(full$secondary_dose))
full$primary_dose <- as.numeric(as.character(full$primary_dose))
full$sec_inf <- as.factor(full$sec_inf)
full$sex <- as.factor(full$sex)
full$prim_inf <- as.factor(full$prim_inf)
full$population <- as.factor(full$population)
full <- full %>% dplyr::mutate(presid = replace_na(presid, 0))
full$prim_inf <- relevel(full$prim_inf, ref = 2)
sinfected <- full %>% filter(sec_inf=="1") 
sinfected$secondary_dose <- as.numeric(as.character(sinfected$secondary_dose))
sinfected$sex <- as.factor(sinfected$sex)
sinfected$prim_inf <- as.factor(sinfected$prim_inf)
sinfected$population <- as.factor(sinfected$population)

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

secresist.u<- lm(sec_resist10~prim_inf:presid, data=sinfected)

# model list
secresist_list <- list(secresist.null=secresist.null, secresist.a=secresist.a, secresist.b=secresist.b, secresist.c=secresist.c, 
                       secresist.d=secresist.d, secresist.e=secresist.e, secresist.f=secresist.f,
                       secresist.g=secresist.g, secresist.h=secresist.h, secresist.i=secresist.i, secresist.j=secresist.j,
                       secresist.k=secresist.k,
                       secresist.l=secresist.l, secresist.m=secresist.m, secresist.n=secresist.n, secresist.o=secresist.o,
                       secresist.p=secresist.p, secresist.q=secresist.q, secresist.r=secresist.r, secresist.s=secresist.s,
                       secresist.t=secresist.t, secresist.u=secresist.u)

# model comparison
aictab(cand.set=secresist_list) 
summary(secresist.j)
#plot(secresist.j)

summary(secresist.null)

# figure
group.colors <- c( "#E69F00","#009E73")
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

# Figure 4
ggplot(sinfected, aes(x = secondary_dose, y = sec_resist10, col=prim_inf)) +
  geom_jitter(width=0.05, height=0, size=2, aes(shape=factor(population))) + 
  scale_x_continuous(trans="log10", breaks=c(30,100, 300,7000))+
  scale_y_continuous(breaks=c(0,1,3,5), labels=c("0", "10", "1000", "100000")) +
  geom_line(data = pred_data, aes(x = secondary_dose, y = pred, col = prim_inf))+
  labs(y= "Second MG load + 1", x = "Second MG dose (CCU/mL)", 
       color="Initial \nsusceptibility") +
  scale_color_manual(values=group.colors, labels=c("Infected", "Uninfected")) +
  theme(legend.position = "right")+
  geom_ribbon(data=pred_data, aes(ymin=lwr.ci, ymax=upr.ci), 
              fill="#000000", alpha=0.2, linetype=0)+
  scale_shape(solid=FALSE)+
  labs(shape="Host \npopulation", colour="Initial \nsusceptibility")


# Second infection resistance (second infected AZ birds)  ------------

sinfected_az <- sinfected %>% filter(population=="AZ") 

# write models
secresist.null2 <- lm(sec_resist10~1, data=sinfected_az)

secresist.a2 <- lm(sec_resist10~secondary_dose, data=sinfected_az) 
secresist.b2 <- lm(sec_resist10~sex, data=sinfected_az)
secresist.c2 <- lm(sec_resist10~prim_inf, data=sinfected_az)
secresist.d2 <- lm(sec_resist10~prim_resist10, data=sinfected_az)
secresist.e2 <- lm(sec_resist10~prim_toleye, data=sinfected_az)
secresist.f2<- lm(sec_resist10~prim_inf:presid, data=sinfected_az)


# model list
secresist_list2 <- list(secresist.null2=secresist.null2, secresist.a2=secresist.a2, 
                       secresist.b2=secresist.b2, secresist.c2=secresist.c2, 
                       secresist.d2=secresist.d2, secresist.e2=secresist.e2, 
                       secresist.f2=secresist.f2)

# model comparison
aictab(cand.set=secresist_list2)

summary(secresist.c2)
summary(secresist.d2)
summary(secresist.e2)
summary(secresist.null2)

#plot(secresist.c2)
#plot(secresist.d2)
#plot(secresist.e2)



# Second infection tolerance (second infected birds) --------

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
sectol.u<- lm(sresid~prim_inf:presid, data=sinfected)

# model list
sectol_list <- list(sectol.null=sectol.null, sectol.a=sectol.a, sectol.b=sectol.b, sectol.c=sectol.c, 
                       sectol.d=sectol.d, sectol.e=sectol.e, sectol.f=sectol.f,
                       sectol.g=sectol.g, sectol.h=sectol.h, sectol.i=sectol.i, sectol.j=sectol.j,
                       sectol.k=sectol.k,
                       sectol.l=sectol.l, sectol.m=sectol.m, sectol.n=sectol.n, sectol.o=sectol.o,
                       sectol.p=sectol.p, sectol.q=sectol.q, sectol.r=sectol.r, sectol.s=sectol.s,
                    sectol.t=sectol.t, sectol.u=sectol.u)

# model comparison
aictab(cand.set=sectol_list)
summary(sectol.null)

summary(sectol.a)
summary(sectol.b)
summary(sectol.u)
summary(sectol.d)

#plot(sectol.a)
#plot(sectol.b)
#plot(sectol.d)
#plot(sectol.u)

# Figure 5
labels<-c("Infected", "Uninfected")

ggplot(sinfected, aes(x = prim_inf, y = sresid)) +
  labs(y= "Second inoculation tolerance", x = "Initial susceptibility", 
       shape="Host \npopulation")+
  theme(legend.position = "right")+
  geom_jitter(width=0.09, height=0, size=3, aes(shape=factor(population)))+
  scale_shape(solid=FALSE) + scale_x_discrete(label=labels) 

sinfected %>% dplyr::group_by(population) %>% summarise(n=n())



# Second infection tolerance (second infected AZ birds) ------------------

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
sectol2.m<- lm(sresid~prim_inf:presid, data=sinfected_az)

# model list
sectol_list2 <- list(sectol2.null=sectol2.null, sectol2.a=sectol2.a,  sectol2.b=sectol2.b, 
                     sectol2.c=sectol2.c, sectol2.d=sectol2.d, sectol2.e=sectol2.e, sectol2.f=sectol2.f,
                     sectol2.g=sectol2.g,sectol2.h=sectol2.h, sectol2.i=sectol2.i, sectol2.j=sectol2.j,
                     sectol2.k=sectol2.k, sectol2.l=sectol2.l, sectol2.m=sectol2.m)

# model comparison
aictab(cand.set=sectol_list2)
summary(sectol2.g)
#plot(sectol2.g)
summary(sectol2.c)
#plot(sectol2.c)
summary(sectol2.null)



