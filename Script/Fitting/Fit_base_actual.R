#Fitting the base model to the actual training data
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#choosing some initial values
initial = c(2.1822324202, #beta0
            0.3886262485, #beta1
            3.9567260886, #phi
            0.0379448249, #nu
            0.0185123191, #p
            0.0069605605, #ra1
            0.0001402397, #ra2
            0.4935457228, #ra3
            0.0003587821, #ra4
            0.1627604763, #ra5
            0.6202412408) #ra6

plot.fit(train.df, initial, known)

#using repeat NM until I get convergence
fit <- hjkb(initial,
            neglk.parameters,
            lower=c(0,0,0,0,0,0,0,0,0,0,0),
            upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
            control=list(info=TRUE),
            data=train.df,
            known=known)

fit.plot <- plot.fit(train.df, fit$par, known)

write.csv(fit$par, 
          file="~/rsv-modelling/Output/Data/fit_base_actual.csv")

capture.output(fit,
               file="~/rsv-modelling/Output/Data/fit_base_actual.txt")

ggsave(filename="Fit_base_actual.png",
       plot=fit.plot,
       path="~/rsv-modelling/Output/Plots",
       width=24,
       height=13.5,
       units="cm",
       dpi=300)

#RSE
RSE <- sum(ssq.parameters(fit$par, data=train.df, known=known)^2) %>% sqrt()
write.csv(RSE, file="~/rsv-modelling/Output/Data/Fit_base_actual_RSE.csv")
