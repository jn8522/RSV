#First round of backwards selection: eliminating one variable
#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")

#reading in the data to fit to
train.df <- read.csv("~/rsv-modelling/Output/Data/train.csv")

#setting these as known values for the SEIR model
known.values = read.csv("~/rsv-modelling/Data/known.values.csv")
known <- known.values$value
names(known) <- known.values$parameter
known <- c(known, new=train.df$cases[1])

#reading in the weather data
weather <- read.csv("~/rsv-modelling/Output/Data/weather.train.csv") %>%
  mutate(precipitation = precipitation/mean(precipitation),
         evaporation = evaporation/mean(evaporation),
         temperature = temperature/mean(temperature),
         dewpoint = dewpoint/mean(dewpoint),
         humidity = humidity/mean(humidity))

#dataframe where we store AIC for each nested model
backwards <- read.csv("~/rsv-modelling/Output/Data/backwards_1.csv")
#computing next fit
include=c(c.p = TRUE,
          c.e = TRUE,
          c.t = FALSE,
          c.d = TRUE,
          c.h = TRUE,
          c.s = TRUE)

#fit to the training data
initial = c(nu = 2.999913e-02,
            p = 1.161561e-02,
            c.0 = 6.466286e-01,
            c.p = 2.918804e-03,
            c.e = -2.738272e-1,
            #c.t = 4.715303e-01,
            c.d = 9.603160e-02,
            c.h = 4.391154e-01,
            c.s = 0.05,
            ra1 = 3.749716e-06,
            ra2 = 1.425406e-05,
            ra3 = 5.927472e-03,
            ra4 = 5.141324e-04,
            ra5 = 3.381306e-01,
            ra6 = 1.010933e-02)

scale = c(1e-7,1e-7,1e-5,
          1e-8, #c.p
          1e-6, #c.e
          #1e-6, #c.t
          1e-6, #c.d
          1e-6, #c.h
          1e-6, #c.s
          1e-11,1e-6,1e-6,1e-9,1e-6,1e-8)

plot.fit.ext(train.df,
             variable=initial,
             known=known,
             weather=weather,
             include=include)

fit <- fit.nm.rpt.ext(initial=initial, 
                      known=known,
                      data=train.df,
                      weather=weather,
                      include=include,
                      print=TRUE,
                      iter=10000,
                      scale=scale)

#calculating AIC
aic <- 14+2*neglk.parameters.ext(fit$par, train.df, known, weather, include)

#storing
backwards[4,3] <- aic

#saving
write.csv(backwards, "~/rsv-modelling/Output/Data/backwards_1.csv")
