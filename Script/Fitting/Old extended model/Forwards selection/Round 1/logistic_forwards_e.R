#Seeing if we can fit the extended model to the training data
library(optimx)
library(dfoptim)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#initial values
initial = c(0.037914307, #nu
            0.018512319, #p
            6, #L
            -1, #c0
            0.25 , #ce
            0.006960561, #ra1
            0.000000000, #ra2
            0.508072090, #ra3
            0.000358782, #ra4
            0.162211160, #ra5
            0.494752960) #ra6

include=c(FALSE, #p
          TRUE, #e
          FALSE, #t
          FALSE, #d
          FALSE, #h
          FALSE) #s

plot.fit.logistic(train.df, initial, known, weather, include)

fit <- hjkb(initial,
            neglk.parameters.logistic,
            lower=c(0,0,0,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

capture.output(fit,
               file="~/rsv-modelling/Output/Data/Forwards_selection/Round 1/logistic_forwards_e.txt")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)

#saving
write.csv(aic, "~/rsv-modelling/Output/Data/Forwards_selection/Round 1/logistic_forwards_e.csv")