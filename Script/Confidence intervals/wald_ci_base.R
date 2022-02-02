#computing Wald CIs for the base model
library(numDeriv)
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#mle for base model
par.base <- read.csv("~/rsv-modelling/Output/Data/base.par.csv", header=FALSE)$V1

#refitting
fit <- hjkb(par.base,
            neglk.parameters,
            lower=c(0,0,0,0,0,0,0,0,0,0,0),
            upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known)

write.csv(fit$par, file="~/rsv-modelling/Output/Data/fit_base_actual.csv")
fit$par-par.base

hessian <- hessian(neglk.parameters,
                   par.base,
                   data=train.df,
                   known=known)
v.hat <- solve(hessian)
wald.intervals <- data.frame(lower=par.base-1.96*sqrt(diag(v.hat)),
                             upper=par.base+1.96*sqrt(diag(v.hat)))