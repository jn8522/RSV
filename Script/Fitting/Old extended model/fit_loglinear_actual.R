#fitting the loglinear model to the training data
#this is the old version that uses six variables
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE_old.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#choosing some initial values by inspection
initial = c(0.0207498229, #nu
            0.0133884709, #p
            0.4586672278, #c.0
            0.0164291838, #c.p
            -0.3074877139, #c.e
            0.7941549076, #c.t
            -0.1650337046, #c.d
            0.6881345737, #c.h
            -0.1210520853, #c.s
            0.0002737532, #ra1
            0.0048107270, #ra2
            0.0206285765, #ra3
            0.0008969523, #ra4
            0.5349117566, #ra5
            0.1446471333) #ra6

include=c(TRUE, #c.p
          TRUE, #c.e
          TRUE, #c.t
          TRUE, #c.d
          TRUE, #c.h
          TRUE) #c.s

fit <- hjkb(initial,
            neglk.parameters.loglinear,
            lower=c(0,0,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0),
            upper=c(Inf,1,Inf,Inf,Inf,Inf,Inf,Inf,Inf,1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known,
            weather=weather,
            include=include)

fit.plot <- plot.fit.loglinear(train.df, fit$par, known, weather, include)

#calculating AIC
aic <- 2*length(initial)+2*neglk.parameters.loglinear(fit$par, train.df, known, weather, include)

ggsave(filename="Fit_extended_actual.png",
       plot=fit.plot,
       path="~/rsv-modelling/Output/Plots",
       width=24,
       height=13.5,
       units="cm",
       dpi=300)

capture.output(fit,
               file="~/rsv-modelling/Output/Data/Fit_extended_actual.txt")
write.csv(aic,
          file="~/rsv-modelling/Output/Data/Fit_extended_actual.csv")


#RSE
RSE <- sum(ssq.parameters.loglinear(fit$par,
                   data=train.df,
                   known=known,
                   weather=weather,
                   include=include)^2) %>% sqrt()
write.csv(RSE, file="~/rsv-modelling/Output/Data/Fit_extended_actual_RSE.csv")