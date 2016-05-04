

## Load Required Packages and Source Functions ##

library(lme4)
library(boot)
source('Model Checks JB.R')
source('Graph GLMMs.R')


## Create Simulated Dataset ##

site <- c('Field','Forest')
plot <- 1:10
rep <- 1:10
trmt <- c('HighN','LowN')

dat <- expand.grid('Site'=site, 'Plot'=plot, 'Rep'=rep,'N'=trmt)
dat$Plot <- paste(dat$Site, dat$Plot, sep=".")
dat$Ht <- rnorm(40, 25, 2) #height is randomly drawn from a normal distr with mean 25
dat[dat$N=='HighN', 'Ht'] <- dat[dat$N=='HighN', 'Ht'] + rnorm(200, 10, 0.5)
#HighN plants are generally taller
dat[dat$Site=='Forest', 'Ht'] <- dat[dat$Site=='Forest', 'Ht'] + rnorm(200, 5, 1)
#Forest plants are generally taller


## Plot the Data ##

plot(dat$N, dat$Ht)
plot(dat$Site, dat$Ht)


## Model Height as a Function of N, Site and Plot ##

mod.ht <- lmer(Ht ~ N + (1|Site) + (1|Plot), data=dat)
#model height as a function of N with site as a random effect


## Model Checks ##

model.check(mod.ht)


## Bootstrapped 95% CI ##

mod.bt <- boot(mod.ht, bootm, R=1000, M=mod.ht)

plot.coef(mod.bt)
