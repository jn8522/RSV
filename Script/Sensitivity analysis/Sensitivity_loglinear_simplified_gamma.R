library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#fitted values
par.loglinear <- read.csv("~/rsv-modelling/Output/Data/Fits/Backwards selection/Loglinear round 1/Fit_loglinear_v_par.csv")$x
include=c(TRUE,TRUE,FALSE,TRUE)

#alpha: 0.5 to 0.8
#sigma: 1/(2/7) to 1/(6/7)
#gamma: 1/(7/7) to 1/(12/7)
#delta: 0.5 to 0.8

#function that fits the base model to training data
#given list containing known values and initial values for fitted parameters
fit.loglinear <- function(known) {
  out <- hjkb(par.loglinear,
              neglk.parameters.loglinear,
              lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
              upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
              control=list(info=TRUE),
              data=train.df,
              known=known,
              weather=weather,
              include=include)
  return(out)
}

#changing values for alpha and re-fitting
list_gamma <- list(c(known[1:5], 1/(12/7), known[7:10]),
                   c(known[1:5], 1/(7/7), known[7:10]))
sensitivity_gamma <- lapply(list_gamma, fit.loglinear)
capture.output(sensitivity_gamma,
               file="~/rsv-modelling/Output/Data/Sensitivity analysis/loglinear.simplified.gamma.txt")
