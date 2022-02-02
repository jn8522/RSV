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

wrapper.hjk <- function(theta,data,known,weather,include) {
  out <- neglk.parameters.loglinear(theta,data,known,weather,include)
    if(is.na(out)){
      return(Inf)
    }
    return(out)
}

fit.hjkb <- hjkb(initial,
                                            wrapper.hjk,
                                            lower=c(0,0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
                                            upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
                                            control=list(info=TRUE),
                                            data=train.df,
                                            known=known,
                                            weather=weather,
                                            include=include)