#benchmarking some different optimisation methods
library(optimx)
library(dfoptim)
library(rbenchmark)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

include=c(TRUE,
          TRUE,
          TRUE,
          FALSE,
          TRUE,
          TRUE)


initial = c(2.999913e-02,
            1.161561e-02,
            6.466286e-01,
            2.918804e-03,
            -2.738272e-1,
            4.715303e-01,
            #c.d = 9.603160e-02,
            4.391154e-01,
            0.05,
            3.749716e-06,
            1.425406e-05,
            5.927472e-03,
            5.141324e-04,
            3.381306e-01,
            1.010933e-02)

wrapper.mads <- function(theta,data,known,weather,include) {
  if(theta[1]<0
     |theta[2]<0|theta[2]>1
     |theta[9]<0|theta[9]>1
     |theta[10]<0|theta[10]>1
     |theta[11]<0|theta[11]>1
     |theta[12]<0|theta[12]>1
     |theta[13]<0|theta[13]>1
     |theta[14]<0|theta[14]>1){
    return(Inf)
  }
  else{
    out <- neglk.parameters.loglinear(theta,data,known,weather,include)
    if(is.na(out)){
      return(Inf)
    }
    return(out)
  }
}
fit.mads <- mads(initial, 
                 wrapper.mads, 
                 control=list(trace=TRUE),
                 data=train.df,
                 known=known,
                 weather=weather,
                 include=include)
