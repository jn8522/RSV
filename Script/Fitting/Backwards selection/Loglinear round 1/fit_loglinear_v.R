#fitting the loglinear model to the training data
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#choosing some initial values by inspection
initial = read.csv("~/rsv-modelling/Output/Data/Fits/Backwards selection/Loglinear round 1/Fit_loglinear_v_par.csv")$x

include=c(TRUE, #c.p
          TRUE, #c.t
          FALSE, #c.v
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
               file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_v_output.txt")
write.csv(aic,
          file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_v_aic.csv")
write.csv(fit$par,
          file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_v_par.csv")


#MSE
mse <- ssq.parameters.loglinear(fit$par,
                   data=train.df,
                   known=known,
                   weather=weather,
                   include=include)
write.csv(mse, file="~/rsv-modelling/Output/Data/Fits/Backwards selection/Fit_loglinear_v_mse.csv")