#First round of backwards selection: eliminating one variable
#benchmarking some different optimisation methods
library(optimx)
library(dfoptim)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#computing next fit
include=c(c.p = TRUE, #p
          c.e = TRUE, #e
          c.t = FALSE, #t
          c.d = TRUE, #d
          c.h = TRUE, #h
          c.s = FALSE) #s

#fit to the training data
initial = c(1.768233e-02, #nu
            1.406493e-02, #p
            1.494280e+00, #c0
            9.912725e-03, #cp
            -2.615582e-01, #ce
            #4.715303e-01, #ct
            2.165783e-01, #cd
            4.059555e-02, #ch
            #0.05, #cs
            9.957309e-04, #ra1
            7.650527e-02, #ra2
            3.678122e-06, #ra3
            6.348221e-04, #ra4
            5.882615e-01, #ra5
            2.678964e-01) #ra6

fit <- hjkb(initial,
            neglk.parameters.loglinear,
            lower=c(0,0,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.loglinear(fit$par, train.df, known, weather, include)

#saving
write.csv(aic, "~/rsv-modelling/Output/Data/Backwards_selection/Round 2/backwards_s_t.csv")
capture.output(fit, file="~/rsv-modelling/Output/Data/Backwards_selection/Round 2/backwards_s_t.txt")
