#First round of backwards selection: eliminating one variable
#benchmarking some different optimisation methods
library(optimx)
library(dfoptim)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#computing next fit
include=c(c.p = FALSE, #p
          c.e = TRUE, #e
          c.t = TRUE, #t
          c.d = TRUE, #d
          c.h = TRUE, #h
          c.s = TRUE) #s

#fit to the training data
initial = c(2.999913e-02, #nu
            1.161561e-02, #p
            6.466286e-01, #c0
            #2.918804e-03, #cp
            -2.738272e-1, #ce
            4.715303e-01, #ct
            9.603160e-02, #cd
            4.391154e-01, #ch
            0.05, #cs
            3.749716e-06, #ra1
            1.425406e-05, #ra2
            5.927472e-03, #ra3
            5.141324e-04, #ra4
            3.381306e-01, #ra5
            1.010933e-02) #ra6

fit <- hjkb(initial,
            neglk.parameters.loglinear,
            lower=c(0,0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.loglinear(fit$par, train.df, known, weather, include)

#saving
write.csv(aic, "~/rsv-modelling/Output/Data/Backwards_selection/backwards_p.csv")
capture.output(fit, file="~/rsv-modelling/Output/Data/Backwards_selection/backwards_p.txt")
