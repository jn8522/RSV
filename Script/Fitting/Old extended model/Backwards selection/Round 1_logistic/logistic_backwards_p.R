#Seeing if we can fit the extended model to the training data
library(optimx)
library(dfoptim)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#initial values
initial = c(1.690834e-02, #nu
            1.426530e-02, #p
            7.638249e+00, #L
            9.350629e-01, #c.0
            #2.430366e-02, #c.p
            -7.549355e-01, #c.e
            5.330510e-01, #c.t
            1.182105e-01, #c.d
            8.810082e-01, #c.h
            -1.532038e+00, #c.s
            6.435438e-06, #ra1
            2.137959e-02, #ra2
            1.155221e-01, #ra3
            6.279353e-04, #ra4
            5.927880e-01, #ra5
            9.216475e-03) #ra6

include=c(FALSE, #p
          TRUE, #e
          TRUE, #t
          TRUE, #d
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
               file="~/rsv-modelling/Output/Data/Backwards_selection/Round 1_logistic/logistic_backwards_p.txt")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)

#saving
write.csv(aic, "~/rsv-modelling/Output/Data/Backwards_selection/Round 1_logistic/logistic_backwards_p.csv")