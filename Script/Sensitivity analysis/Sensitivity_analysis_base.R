library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()
indices <- data.frame(parameter=read.csv("~/rsv-modelling/Data/known.values.csv")$parameter,
                      index=1:9,
                      fixed.value = known[1:9])

#reading in fitted values
par.base <- read.csv("~/rsv-modelling/Output/Data/base.par.csv", header=FALSE)$V1

#alpha: 0.5 to 0.8
#sigma: 1/(2/7) to 1/(6/7)
#gamma: 1/(7/7) to 1/(12/7)
#delta: 0.5 to 0.8

#function that fits the base model to training data
#given list containing known values and initial values for fitted parameters
fit.base <- function(known) {
  out <- hjkb(par.base,
              neglk.parameters,
              lower=c(0,0,0,0,0,0,0,0,0,0,0),
              upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
              control=list(info=TRUE),
              data=train.df,
              known=known)
  return(out)
}

# #changing values for alpha and re-fitting
# list_alpha <- list(c(known[1:3], 0.5, known[5:10]),
#                    c(known[1:3], 0.8, known[5:10]))
# sensitivity_alpha <- lapply(list_alpha, fit.base)
# capture.output(sensitivity_alpha,
#                file="~/rsv-modelling/Output/Data/Sensitivity analysis/base.alpha.txt")
# 
# #changing values for delta and re-fitting
# list_delta <- list(c(known[1:6], 0.5, known[8:10]),
#                    c(known[1:6], 0.8, known[8:10]))
# sensitivity_delta <- lapply(list_delta, fit.base)
# capture.output(sensitivity_delta,
#                file="~/rsv-modelling/Output/Data/Sensitivity analysis/base.delta.txt")

#changing values for gamma and re-fitting
list_gamma <- list(c(known[1:5], 1/(12/7), known[7:10]),
                   c(known[1:5], 1/(7/7), known[7:10]))
sensitivity_gamma <- lapply(list_gamma, fit.base)
capture.output(sensitivity_gamma,
               file="~/rsv-modelling/Output/Data/Sensitivity analysis/base.gamma.txt")

#changing values for sigma and re-fitting
list_sigma <- list(c(known[1:4], 1/(6/7), known[6:10]),
                   c(known[1:4], 1/(2/7), known[6:10]))
sensitivity_sigma <- lapply(list_sigma, fit.base)
capture.output(sensitivity_sigma,
               file="~/rsv-modelling/Output/Data/Sensitivity analysis/base.sigma.txt")
