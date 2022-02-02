#Seeing if we can fit the logistic extended model to the training data
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#choosing some initial values by inspection
initial = c(0.010089905, #nu
            0.023246359, #p
            6.213323199, #L
            -4.736798317, #c.0
            1.666597761, #c.p
            #, #c.t
            6.211416113, #c.v
            -5.880133346, #c.s
            0.000727346, #ra1
            0.581840515, #ra2
            6.86646E-05, #ra3
            0.00000803, #ra4
            0.603408814, #ra5
            0.021545411) #ra6

include=c(TRUE, #c.p
          FALSE, #c.t
          TRUE, #c.v
          TRUE) #c.s

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
               file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_logistic_t_output.txt")
write.csv(fit$par,
          file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_logistic_t_par.csv")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)
write.csv(aic, "~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_logistic_t_aic.csv")

#mse
mse <- ssq.parameters.logistic(fit$par,
                                    data=train.df,
                                    known=known,
                                    weather=weather,
                                    include=include)
write.csv(mse, file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_logistic_t_mse.csv")