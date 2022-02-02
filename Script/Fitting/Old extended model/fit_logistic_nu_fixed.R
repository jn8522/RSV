#seeing how the logistic model's fit changes if we fix nu to give a short period of immunity
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE_old.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#re-defining these functions from ODE.R so that nu is fixed
neglk.parameters.logistic.fixed <- function(theta,
                                      data, #observed df
                                      known,
                                      weather,
                                      include){
  #solving ODE for the given parameters
  fitted.df <- ode.solve.logistic(variable=c(0.0379143073, theta), 
                                  known=known,
                                  noise=FALSE,
                                  weather=weather,
                                  include)%>%
    filter(week>1)
  
  observed.df <- data %>% filter(week>1)
  #calculating likelihood of observed data
  p = theta[1]
  lk <- -sum(dnorm(observed.df$cases,
                   mean = fitted.df$cases,
                   sd = sqrt(fitted.df$cases*(1-p)),
                   log = TRUE))
  if(is.na(lk)){
    return(Inf)
  }
  else{
    return(lk)
  }
}

initial <- c(0.01379609, #p
             6.28315, #L
             -48.387, #c.0
             0.2102167, #c.p
             35.53598, #c.t
             -20.18632, #c.d
             35.72293, #c.h
             -1.519312, #c.s
             0.003630398, #ra1
             0.2559682, #ra2
             9.15527E-05, #ra3
             7.48042E-05, #ra4
             0.461448, #ra5
             0.009155273) #ra6
include=c(TRUE,FALSE,TRUE,TRUE,TRUE,TRUE)

fit <- hjkb(initial,
            neglk.parameters.logistic.fixed,
            lower=c(0,0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(1,Inf,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

capture.output(fit,
               file="~/rsv-modelling/Output/Data/Fits/logistic_nu_fixed.txt")