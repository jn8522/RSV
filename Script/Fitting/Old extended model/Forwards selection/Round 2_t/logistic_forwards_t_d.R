#Seeing if we can fit the extended model to the training data
library(optimx)
library(dfoptim)

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#initial values
initial = c(2.007678e-02, #nu
            1.340825e-02, #p
            8.076458e+00, #L
            1.358200e+00, #c0
            -1.156876e+00, #ct
            0, #cd
            7.690430e-03, #ra1
            4.724121e-02,
            2.054749e-01,
            0.000000e+00,
            5.264435e-01,
            1.525879e-05) #ra6

include=c(FALSE, #p
          FALSE, #e
          TRUE, #t
          TRUE, #d
          FALSE, #h
          FALSE) #s

plot.fit.logistic(train.df, initial, known, weather, include)

fit <- hjkb(initial,
            neglk.parameters.logistic,
            lower=c(0,0,0,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

capture.output(fit,
               file="~/rsv-modelling/Output/Data/Forwards_selection/Round 2_t/logistic_forwards_t_d.txt")

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.logistic(fit$par, train.df, known, weather, include)

#saving
write.csv(aic, "~/rsv-modelling/Output/Data/Forwards_selection/Round 2_t/logistic_forwards_t_d.csv")