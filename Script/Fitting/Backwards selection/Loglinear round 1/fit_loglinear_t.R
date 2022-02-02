#fitting the loglinear model to the training data
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#choosing some initial values by inspection
initial = c(0.033824952, #nu
            0.011035512, #p
            1.290411585, #c.0
            0.064069429, #c.p
            #-0.499265105, #c.t
            0.024017361, #c.v
            0.045433055, #c.s
            0.007815053, #ra1
            0.000183105, #ra2
            0.106933594, #ra3
            0, #ra4
            0.28478986, #ra5
            0.644105987) #ra6

include=c(TRUE, #c.p
          FALSE, #c.t
          TRUE, #c.v
          TRUE) #c.s

plot.fit.loglinear(train.df,
                   initial,
                   known,
                   weather,
                   include)

fit <- hjkb(initial,
            neglk.parameters.loglinear,
            lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.loglinear(fit$par, train.df, known, weather, include)

capture.output(fit,
               file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_t_output.txt")
write.csv(aic,
          file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_t_aic.csv")
write.csv(fit$par,
          file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_t_par.csv")


#MSE
mse <- ssq.parameters.loglinear(fit$par,
                   data=train.df,
                   known=known,
                   weather=weather,
                   include=include)
write.csv(mse, file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_t_mse.csv")