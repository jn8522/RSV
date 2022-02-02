library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#fitted values
par.logistic <- read.csv("~/rsv-modelling/Output/Data/extended.par.csv")$x
include=c(TRUE,FALSE,TRUE,TRUE,TRUE,TRUE)

#alpha: 0.5 to 0.8
#sigma: 1/(2/7) to 1/(6/7)
#gamma: 1/(7/7) to 1/(12/7)
#delta: 0.5 to 0.8

#function that fits the base model to training data
#given list containing known values and initial values for fitted parameters
fit.log <- function(known) {
  out <- hjkb(par.logistic,
              neglk.parameters.logistic,
              lower=c(0,0,0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
              upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
              control=list(info=TRUE),
              data=train.df,
              known=known,
              weather=weather,
              include=include)
  return(out)
}

#changing values for alpha and re-fitting
list_alpha <- list(c(known[1:3], 0.5, known[5:10]),
                   c(known[1:3], 0.8, known[5:10]))
sensitivity_alpha <- lapply(list_alpha, fit.log)
capture.output(sensitivity_alpha,
               file="~/rsv-modelling/Output/Data/Sensitivity analysis/logistic.alpha.txt")
