#Seeing if we can fit the extended model to the training data
library(optimx)
library(dfoptim)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#initial values
initial = c(1.686638e-02, #nu
            1.379609e-02, #p
            6.283150e+00, #L
            0, #c0
            2.102167e-01, #cp
            3.553598e+01, #ct
            -2.018632e+01, #cd
            #3.572293e+01, #ch
            -1.519312e+00, #cs
            3.630398e-03, #ra1
            2.559682e-01, #ra2
            9.155273e-05, #ra3
            7.480420e-05, #ra4
            4.614480e-01, #ra5
            9.155273e-03) #ra6

include=c(TRUE, #p
          FALSE, #e
          TRUE, #t
          TRUE, #d
          FALSE, #h
          TRUE) #s

#using repeat nm
fit <- hjkb(initial,
            neglk.parameters.logistic,
            lower=c(0,0,0,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

capture.output(fit,
               file="~/rsv-modelling/Output/Data/Backwards_selection/Round 2_logistic/logistic_backwards_e_h.txt")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)

#saving
write.csv(aic, "~/rsv-modelling/Output/Data/Backwards_selection/Round 2_logistic/logistic_backwards_e_h.csv")