#Seeing if we can fit the extended model to the training data
library(optimx)
library(dfoptim)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#initial values
initial = c(0.0231108693, #nu
            0.0124699184, #p
            5.9680164510, #L
            0.1876741723, #c.0
            0.6423432668, #c.p
            -0.4575480227, #c.e
            0.7765060298, #c.t
            #1.182105e-01, #c.d
            0.2518679071, #c.h
            -0.5819287278, #c.s
            0.0012576832, #ra1
            0.0069250578, #ra2
            0.0024062569, #ra3
            0.0004271513, #ra4
            0.5086949769, #ra5
            0.9999483309) #ra6

include=c(TRUE, #p
          TRUE, #e
          TRUE, #t
          FALSE, #d
          TRUE, #h
          TRUE) #s

#using repeat nm
fit <- hjkb(initial,
            neglk.parameters.logistic,
            lower=c(0,0,0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

capture.output(fit,
               file="~/rsv-modelling/Output/Data/Backwards_selection/Round 1_logistic/logistic_backwards_d.txt")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)

#saving
write.csv(aic, "~/rsv-modelling/Output/Data/Backwards_selection/Round 1_logistic/logistic_backwards_d.csv")