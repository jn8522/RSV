#Seeing if we can fit the logistic extended model to the training data
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#choosing some initial values by inspection
initial = c(0.011005432, #nu
            0.022800039, #p
            5.831708514, #L
            -4.789593727, #c.0
            1.717745222, #c.p
            -0.493169936, #c.t
            6.262197363, #c.v
            #-5.670935348, #c.s
            0.001101186, #ra1
            0.656272888, #ra2
            3.05E-05, #ra3
            8.03E-06, #ra4
            0.582427979, #ra5
            0.000572205) #ra6

include=c(TRUE, #c.p
          TRUE, #c.t
          TRUE, #c.v
          FALSE) #c.s

fit <- hjkb(initial,
            neglk.parameters.logistic,
            lower=c(0,0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

capture.output(fit,
               file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_logistic_s_output.txt")
write.csv(fit$par,
          file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_logistic_s_par.csv")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)
write.csv(aic, "~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_logistic_s_aic.csv")

#mse
mse <- ssq.parameters.logistic(fit$par,
                                    data=train.df,
                                    known=known,
                                    weather=weather,
                                    include=include)
write.csv(mse, file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_logistic_s_mse.csv")