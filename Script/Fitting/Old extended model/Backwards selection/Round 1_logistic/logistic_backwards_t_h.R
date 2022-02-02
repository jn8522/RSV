#Seeing if we can fit the extended model to the training data
library(optimx)
library(dfoptim)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#initial values
initial = c(2.029835e-02, #nu
            1.360014e-02, #p
            9.492116e+00, #L
            1.191765e-01, #c0
            3.002137e-03, #cp
            -5.393523e-01, #ce
            3.073192e-01, #cd
            -1.768819e-01, #cs
            1.663766e-05, #ra1
            3.730597e-01,
            4.234464e-01,
            1.464227e-03,
            5.532231e-01,
            5.825607e-05) #ra6

include=c(TRUE, #p
          TRUE, #e
          FALSE, #t
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
               file="~/rsv-modelling/Output/Data/Backwards_selection/Round 1_logistic/logistic_backwards_t_h.txt")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)

#saving
write.csv(aic, "~/rsv-modelling/Output/Data/Backwards_selection/Round 1_logistic/logistic_backwards_t_h.csv")