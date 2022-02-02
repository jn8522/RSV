#fitting the loglinear model to the training data
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#choosing some initial values by inspection
initial = c(0.034305604, #nu
            0.010924886, #p
            1.361334437, #c.0
            0.081525484, #c.p
            -0.327878386, #c.t
            0, #c.v
            0.101341258, #c.s
            0.009284973, #ra1
            0, #ra2
            0.275390625, #ra3
            0, #ra4
            0.307258606, #ra5
            0.048339844) #ra6

include=c(TRUE, #c.p
          TRUE, #c.t
          TRUE, #c.v
          TRUE) #c.s

plot.fit.loglinear(train.df,
                   initial,
                   known,
                   weather,
                   include)

# fit <- hjkb(initial,
#             neglk.parameters.loglinear,
#             lower=c(0,0,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
#             upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
#             control=list(info=TRUE),
#             data=train.df,
#             known=known,
#             weather=weather,
#             include=include)

wrapper.mads <- function(theta,data,known,weather,include) {
        if(theta[1]<0
           |theta[2]<0|theta[2]>1
           |theta[8]<0|theta[8]>1
           |theta[9]<0|theta[9]>1
           |theta[10]<0|theta[10]>1
           |theta[11]<0|theta[11]>1
           |theta[12]<0|theta[12]>1
           |theta[13]<0|theta[13]>1){
                return(Inf)
        }
        else{
                out <- neglk.parameters.loglinear(theta,data,known,weather,include)
                if(is.na(out)){
                        return(Inf)
                }
                return(out)
        }
}
fit <- mads(initial, 
            wrapper.mads, 
            control=list(trace=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)


#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.loglinear(fit$par, train.df, known, weather, include)


capture.output(fit,
               file="~/rsv-modelling/Output/Data/Fits/Fit_loglinear_output.txt")
write.csv(aic,
          file="~/rsv-modelling/Output/Data/Fits/Fit_loglinear_aic.csv")
write.csv(fit$par,
          file="~/rsv-modelling/Output/Data/Fits/Fit_loglinear_par.csv")


#RSE
mse <- ssq.parameters.loglinear(fit$par,
                   data=train.df,
                   known=known,
                   weather=weather,
                   include=include)
write.csv(mse, file="~/rsv-modelling/Output/Data/Fits/Fit_loglinear_mse.csv")