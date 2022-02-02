#fitting the loglinear model to the training data
library(numDeriv)
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#mle for simplified loglinear
par.loglinear <- read.csv("~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_v_par.csv")$x
include=c(TRUE,TRUE,FALSE,TRUE)

hessian <- hessian(neglk.parameters.loglinear,
                       par.loglinear,
                       data=train.df,
                       known=known,
                       weather=weather,
                       include=include)
v.hat <- solve(hessian)
wald.intervals <- data.frame(lower=par.loglinear-1.96*sqrt(diag(v.hat)),
                             upper=par.loglinear+1.96*sqrt(diag(v.hat)))