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
          c.t = TRUE, #t
          c.d = TRUE, #d
          c.h = FALSE, #h
          c.s = FALSE) #s

#fit to the training data
initial = c(0.0224585753, #nu
            0.0129782421, #p
            1.3788822747, #c0
            0.0150696044, #cp
            -0.3166217232, #ce
            0.0618367536, #ct
            0.2168520144, #cd
            #4.391154e-01, #ch
            #0.05, #cs
            0.0006597589, #ra1
            0.0652784542, #ra2
            0.0007780459, #ra3
            0.0005974058, #ra4
            0.4939003416, #ra5
            0.2202802728) #ra6

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
write.csv(aic, "~/rsv-modelling/Output/Data/Backwards_selection/Round 2/backwards_s_h.csv")
capture.output(fit, file="~/rsv-modelling/Output/Data/Backwards_selection/Round 2/backwards_s_h.txt")
