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
          c.h = FALSE, #h
          c.s = TRUE) #s

#fit to the training data
initial = c(2.601099e-02, #nu
            1.231825e-02, #p
            1.487233e+00, #c0
            2.406154e-03, #cp
            -3.596435e-01, #ce
            #4.715303e-01, #ct
            2.713719e-01, #cd
            #4.391154e-01, #ch
            -1.671310e-01, #cs
            3.569718e-06, #ra1
            5.026222e-02, #ra2
            5.802667e-03, #ra3
            7.173295e-04, #ra4
            4.272169e-01, #ra5
            9.123048e-03) #ra6

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
write.csv(aic, "~/rsv-modelling/Output/Data/Backwards_selection/Round 2/backwards_t_h.csv")
capture.output(fit, file="~/rsv-modelling/Output/Data/Backwards_selection/Round 2/backwards_t_h.txt")
