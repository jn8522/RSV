#Here we generate some data from the base SEIR model and see if we can fit to it
source("~/rsv-modelling/Script/Functions/ODE.R")
library(minpack.lm)
set.seed(1)

#setting these as known values for the SEIR model
known.values = read.csv("~/rsv-modelling/Data/known.values.csv")
known <- known.values$value
names(known) <- known.values$parameter

#generating some data with these parameters and initial conditions
variable.parms.generate = c(beta0 = 1.99,
                            beta1 = 0.65,
                            phi = 2.43,
                            nu = 1/23.5,
                            p=0.05,
                            ra1=0.5,
                            ra2=0.5,
                            ra3=0.5,
                            ra4=0.5,
                            ra5=0.5,
                            ra6=0.5)

generated.df <- ode.solve(variable=variable.parms.generate, 
                          known=c(known, new=0),
                          weeks=15*52,
                          noise=TRUE) %>%
  filter(week>(5*52)) %>%
  mutate(week=week-5*52)

ggplot(generated.df, aes(x=week, y=cases))+geom_line()

#now we try fitting to this fake data
#update known parameters to include number of new cases in the first week
known = c(known, new=generated.df$cases[1])

#choose something vaguely sensible by inspecting plot
initial = c(beta0 = 1.9,
            beta1 = 0.6,
            phi = 2.4,
            nu = 1/20,
            p=0.07,
            ra1=0.01,
            ra2=0.3,
            ra3=0.6,
            ra4=0.01,
            ra5=0.3,
            ra6=0.6)
plot.fit(generated.df,
         variable=initial,
         known=known)

#fitting with nelder mead
fit.nm <- optim(par=initial,
      fn = neglk.parameters,
      control = list(trace=TRUE,
                     maxit=10000,
                     parscale=c(1/170143,
                                  1/392045,
                                  1/127266.4,
                                  1/7620024,
                                  1/4421041,
                                  1/29719720,
                                  1/752537.7,
                                  1/356297.4,
                                  1/30055285,
                                  1/1329952,
                                  1/465658.4)),
      known=known,
      data=generated.df)
fit.nm
plot.fit(generated.df,
         variable=fit.nm$par,
         known=known)

#fitting with LM
fit.nls <- nls.lm(fit.nm$par,
                  lower=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                  upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1),
                  fn = ssq.parameters,
                  data=generated.df,
                  known=known)
summary(fit.nls)
plot.fit(generated.df,
         variable=fit.nls$par,
         known=known)