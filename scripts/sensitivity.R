library(dfoptim)
library(parallel)
source("src/fit.R")
source("src/ODE.R")
source("src/plot.R")

#reading in the model fits
fit.base <- readRDS("models/fit.base.Rdata")
fit.weather <- readRDS("models/fit.weather.Rdata")

#reading in some other relevant data
RSV <- read.csv("data/processed/RSV.csv")
known <- c(read.csv("data/raw/known.values.csv")$value, RSV$cases[1])
weather.school <- read.csv("data/processed/weather.school.csv")%>%filter(week<=nrow(RSV))
include=c(TRUE, #c.p
          TRUE, #c.t
          FALSE, #c.v
          TRUE) #c.s

#making a list of all the vectors of known parameters we want to fit over
known.matrix <- matrix(rep(known,5),nrow=10,ncol=5,byrow=FALSE)
rownames(known.matrix) <- c(read.csv("data/raw/known.values.csv")$parameter, "y0")

alpha.matrix <- known.matrix
alpha.matrix["alpha",]<-seq(0.5,0.8,length.out=5)

delta.matrix <- known.matrix
delta.matrix["delta",]<-seq(0.5,0.8,length.out=5)

gamma.matrix <- known.matrix
gamma.matrix["gamma",]<-seq(1/(12/7),1/(7/7),length.out=5)

sigma.matrix <- known.matrix
sigma.matrix["sigma",]<-seq(1/(6/7),1/(2/7),length.out=5)

knowns.list <- list()
for (i in 1:5) {
  knowns.list[[i]] <- alpha.matrix[,i]
  knowns.list[[i+5]] <- sigma.matrix[,i]
  knowns.list[[i+10]] <- gamma.matrix[,i]
  knowns.list[[i+15]] <- delta.matrix[,i]
}
#BASE MODEL
#fitting over them
fits.sensitivity.base <- mcmapply(hjkb,
                                  mc.cores = 11,
                                  known=knowns.list,
                                  MoreArgs=list(fn=neglk.parameters.base,
                                                lower=c(0,0,0,0,0,0,0,0,0,0,0),
                                                upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                                                control=list(info=TRUE),
                                                par=fit.base$par,
                                                data=RSV,
                                                i=0,
                                                k=0))

#checking that they really converged
fits.sensitivity.base2 <- mcmapply(hjkb,
                                  mc.cores = 11,
                                  par=fits.sensitivity.base[1,],
                                  known=knowns.list,
                                  MoreArgs=list(fn=neglk.parameters.base,
                                                lower=c(0,0,0,0,0,0,0,0,0,0,0),
                                                upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                                                control=list(info=TRUE),
                                                data=RSV,
                                                i=0,
                                                k=0))
saveRDS(fits.sensitivity.base2, file="models/fit.sensitivity.base.RData")

#WEATHER MODEL
#fitting over them
fits.sensitivity.weather <- mcmapply(hjkb,
                                  mc.cores = 11,
                                  known=knowns.list,
                                  MoreArgs=list(fn=neglk.parameters.weather,
                                                lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                                                upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                                                control=list(info=TRUE),
                                                par=fit.weather$par,
                                                data=RSV,
                                                i=0,
                                                k=0,
                                                weather=weather.school,
                                                include=include))

#checking that they really converged
fits.sensitivity.weather2 <- mcmapply(hjkb,
                                   mc.cores = 11,
                                   par=fits.sensitivity.weather[1,],
                                   known=knowns.list,
                                   MoreArgs=list(fn=neglk.parameters.weather,
                                                 lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                                                 upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                                                 control=list(info=TRUE),
                                                 data=RSV,
                                                 i=0,
                                                 k=0,
                                                 weather=weather.school,
                                                 include=include))
saveRDS(fits.sensitivity.weather2, file="models/fits.sensitivity.weather.Rdata")