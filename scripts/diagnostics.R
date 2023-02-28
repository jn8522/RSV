#checking how sensitive the fit of the base model is to choice of starting parameters
library(dfoptim)
library(parallel)

source("src/fit.R")
source("src/ODE.R")

#reading in data
RSV <- read.csv("data/processed/RSV.csv")
known <- c(read.csv("data/raw/known.values.csv")$value, RSV$cases[1])
names(known) <- c("mu","eta1","eta2","alpha","sigma","gamma","delta","n1","n2","y0")
weather.school <- read.csv("data/processed/weather.school.csv")%>%filter(week<=nrow(RSV))
include=c(TRUE, #c.p
          TRUE, #c.t
          FALSE, #c.v
          TRUE) #c.s

#randomly generating some starting parameters
set.seed(1)
initials <- list()
for (i in 1:4) {
  initials[[i]] = c(beta0=runif(1,0,10),
                    beta1=runif(1,0,1),
                    phi=runif(1,0,2*pi),
                    nu=runif(1,0.1,1),
                    p=runif(1,0,1),
                    ra1=runif(1,0,1),
                    ra2=runif(1,0,1),
                    ra3=runif(1,0,1),
                    ra4=runif(1,0,1),
                    ra5=runif(1,0,1),
                    ra6=runif(1,0,1))
}

diagnostic.fits.base <- mcmapply(hjkb,
                                 mc.cores = 4,
                                 par=initials,
                                 MoreArgs=list(fn=neglk.parameters.base,
                                               lower=c(0,0,0,0,0,0,0,0,0,0,0),
                                               upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                                               control=list(info=TRUE),
                                               data=RSV,
                                               known=known,
                                               i=0, #fixing i=k=0 so we're not using cross-validation
                                               k=0))
#checking that they definitely converged
diagnostic.fits.base2 <- mcmapply(hjkb,
                                   mc.cores = 4,
                                   par=diagnostic.fits.base["par",],
                                   MoreArgs=list(fn=neglk.parameters.base,
                                                 lower=c(0,0,0,0,0,0,0,0,0,0,0),
                                                 upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                                                 control=list(info=TRUE),
                                                 data=RSV,
                                                 known=known,
                                                 i=0, #fixing i=k=0 so we're not using cross-validation
                                                 k=0))

#storing the result
saveRDS(diagnostic.fits.base2, file="models/diagnostic.fits.base.Rdata")

#and for the weather model
set.seed(2)
initials.weather <- list()
for (i in 1:4) {
  initials.weather[[i]] = c(nu=runif(1,0.01,1), #nu
                            p=runif(1,0,1), #p
                            c.0=runif(1,1,2), #c.0
                            c.p=runif(1,0,1), #c.p
                            c.t=runif(1,-1,0), #c.t
                            c.s=runif(1,0,5), #c.s
                            ra1=runif(1,0.1,1), #ra1
                            ra2=runif(1,0,1), #ra2
                            ra3=runif(1,0,1), #ra3
                            ra4=runif(1,0,1), #ra4
                            ra5=runif(1,0,1), #ra5
                            ra6=runif(1,0,1)) #ra6
}

diagnostic.fits.weather <- mcmapply(hjkb,
                                    mc.cores = 4,
                                    par=initials.weather,
                                    MoreArgs=list(fn=neglk.parameters.weather,
                                                  lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                                                  upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                                                  control=list(info=TRUE),
                                                  data=RSV,
                                                  known=known,
                                                  weather=weather.school,
                                                  include=include,
                                                  i=0, #fixing i=k=0 so we're not using cross-validation
                                                  k=0))
#checking that they definitely converged
diagnostic.fits.weather2 <- mcmapply(hjkb,
                                    mc.cores = 4,
                                    par=diagnostic.fits.weather["par",],
                                    MoreArgs=list(fn=neglk.parameters.weather,
                                                  lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                                                  upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                                                  control=list(info=TRUE),
                                                  data=RSV,
                                                  known=known,
                                                  weather=weather.school,
                                                  include=include,
                                                  i=0, #fixing i=k=0 so we're not using cross-validation
                                                  k=0))

#storing the result
saveRDS(diagnostic.fits.weather2, file="models/diagnostic.fits.weather.Rdata")