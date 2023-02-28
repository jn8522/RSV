library(dfoptim)
library(parallel)
source("src/fit.R")
source("src/ODE.R")
source("src/plot.R")

#reading in data
RSV <- read.csv("data/processed/RSV.csv")
known <- c(read.csv("data/raw/known.values.csv")$value, RSV$cases[1])
weather.school <- read.csv("data/processed/weather.school.csv")%>%filter(week<=nrow(RSV))
include=c(TRUE, #c.p
          TRUE, #c.t
          FALSE, #c.v
          TRUE) #c.s

#BASE MODEL
#choosing some initial values
initial = c(2.1822324202, #beta0
            0.3886262485, #beta1
            3.9567260886, #phi
            0.0379448249, #nu
            0.0185123191, #p
            0.0069605605, #ra1
            0.0001402397, #ra2
            0.4935457228, #ra3
            0.0003587821, #ra4
            0.1627604763, #ra5
            0.6202412408) #ra6

#fitting the base model using 5-fold cross-validation
X <- c(i=1:5)
fits.base.CV <- mclapply(X, hjkb,
                      mc.cores=5,
                      par=initial,
                      fn=neglk.parameters.base,
                      lower=c(0,0,0,0,0,0,0,0,0,0,0),
                      upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                      control=list(info=TRUE),
                      data=RSV,
                      known=known,
                      k=5)

#fitting one more time to check that it *really* converged
fits.base.CV2 <- mcmapply(hjkb,
                          mc.cores = 5,
                          par=list(fits.base.CV$i1$par,
                                   fits.base.CV$i2$par,
                                   fits.base.CV$i3$par,
                                   fits.base.CV$i4$par,
                                   fits.base.CV$i5$par),
                          i=1:5,
                          MoreArgs=list(fn=neglk.parameters.base,
                                        lower=c(0,0,0,0,0,0,0,0,0,0,0),
                                        upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                                        control=list(info=TRUE),
                                        data=RSV,
                                        known=known,
                                        k=5))

saveRDS(fits.base.CV2, file="models/fits.base.CV.Rdata")

#fitting the base model to all of the data (training and test) -- this is needed for figures 3-5
fit.base <- hjkb(fits.base.CV2[1,1]$par,
                 neglk.parameters.base,
                 lower=c(0,0,0,0,0,0,0,0,0,0,0),
                 upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                 control=list(info=TRUE),
                 data=RSV,
                 known=known,
                 i=0,
                 k=0)
#checking that it definitely converged
fit.base2 <- hjkb(fit.base$par,
                 neglk.parameters.base,
                 lower=c(0,0,0,0,0,0,0,0,0,0,0),
                 upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                 control=list(info=TRUE),
                 data=RSV,
                 known=known,
                 i=0,
                 k=0)
names(fit.base2$par)<-c("beta0",
                       "beta1",
                       "phi",
                       "nu",
                       "p",
                       "ra1",
                       "ra2",
                       "ra3",
                       "ra4",
                       "ra5",
                       "ra6")
saveRDS(fit.base2, file="models/fit.base.Rdata")

#WEATHER MODEL
#choosing some initial values
initial.weather = c(0.03, #nu
            0.01, #p
            1.36, #c.0
            0.08, #c.p
            -0.32, #c.t
            0.10, #c.s
            0.009, #ra1
            0, #ra2
            0.27, #ra3
            0, #ra4
            0.30, #ra5
            0.04) #ra6

X <- c(i=1:5)
fits.weather.CV <- mclapply(X, hjkb,
                            mc.cores=5,
                            par=initial.weather,
                            fn=neglk.parameters.weather,
                            lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                            upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                            control=list(info=TRUE),
                            data=RSV,
                            known=known,
                            weather=weather.school,
                            include=include,
                            k=5)

#fitting one more time to check that it *really* converged
fits.weather.CV2 <- mcmapply(hjkb,
                             mc.cores = 5,
                             par=list(fits.weather.CV$i1$par,
                                      fits.weather.CV$i2$par,
                                      fits.weather.CV$i3$par,
                                      fits.weather.CV$i4$par,
                                      fits.weather.CV$i5$par),
                             i=1:5,
                             MoreArgs=list(fn=neglk.parameters.weather,
                                           lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                                           upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                                           control=list(info=TRUE),
                                           data=RSV,
                                           known=known,
                                           weather=weather.school,
                                           include=include,
                                           k=5))

saveRDS(fits.weather.CV2, file="models/fits.weather.CV.Rdata")

#fitting the weather model to all of the data (training and test pooled) -- this is needed to produce figures 3-5
fit.weather <- hjkb(fits.weather.CV2[1,1]$par,
                    neglk.parameters.weather,
                    lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                    upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                    control=list(info=TRUE),
                    data=RSV,
                    known=known,
                    weather=weather.school,
                    include=include,
                    i=0,
                    k=0)
#checking that it definitely converged
fit.weather2 <- hjkb(fit.weather$par,
                     neglk.parameters.weather,
                     lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                     upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                     control=list(info=TRUE),
                     data=RSV,
                     known=known,
                     weather=weather.school,
                     include=include,
                     i=0,
                     k=0)
names(fit.weather2$par) <- c("nu",
                             "p",
                             "c.0",
                             "c.p",
                             "c.t",
                             "c.s",
                             "ra1",
                             "ra2",
                             "ra3",
                             "ra4",
                             "ra5",
                             "ra6")
saveRDS(fit.weather2, file="models/fit.weather.Rdata")
