#Seeing if we can fit the logistic extended model to the training data
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#choosing some initial values by inspection
initial = c(0.011005432, #nu
            0.016128133, #p
            9.149091326, #L
            4.616131813, #c.0
            0.961885847, #c.p
            -0.49730858, #c.t
            0, #c.v
            -4.036902145, #c.s
            0.008913686, #ra1
            0.39074707, #ra2
            0.661376953, #ra3
            8.03E-06, #ra4
            0.725860596, #ra5
            0.905029297) #ra6

include=c(TRUE, #c.p
          TRUE, #c.t
          TRUE, #c.v
          TRUE) #c.s

plot.fit.logistic(train.df,
         variable=initial,
         known=known,
         weather=weather,
         include=include)

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
               file="~/rsv-modelling/Output/Data/Fits/Fit_logistic_output.txt")
write.csv(fit$par,
          file="~/rsv-modelling/Output/Data/Fits/Fit_logistic_par.csv")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)
write.csv(aic, "~/rsv-modelling/Output/Data/Fits/Fit_logistic_aic.csv")

#mse
mse <- ssq.parameters.logistic(fit$par,
                                    data=train.df,
                                    known=known,
                                    weather=weather,
                                    include=include)
write.csv(mse, file="~/rsv-modelling/Output/Data/Fits/Fit_logistic_mse.csv")