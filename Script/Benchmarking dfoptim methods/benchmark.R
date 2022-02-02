#seeing what the slow steps are in ode.solve.loglinear
library(profvis)
source("~/rsv-modelling/Script/Functions/ODE.R")

#usual preamble
#reading in training data
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

include=c(c.p = TRUE,
          c.e = TRUE,
          c.t = TRUE,
          c.d = FALSE,
          c.h = TRUE,
          c.s = FALSE)
initial = c(nu = 2.999913e-02,
            p = 1.161561e-02,
            c.0 = 6.466286e-01,
            c.p = 2.918804e-03,
            c.e = -2.738272e-1,
            c.t = 4.715303e-01,
            #c.d = 9.603160e-02,
            c.h = 4.391154e-01,
            #c.s = 0.05,
            ra1 = 3.749716e-06,
            ra2 = 1.425406e-05,
            ra3 = 5.927472e-03,
            ra4 = 5.141324e-04,
            ra5 = 3.381306e-01,
            ra6 = 1.010933e-02)
scale = c(1e-7,1e-7,1e-5,
          1e-8, #c.p
          1e-6, #c.e
          1e-6, #c.t
          #1e-6, #c.d
          1e-6, #c.h
          #1e-6, #c.s
          1e-11,1e-6,1e-6,1e-9,1e-6,1e-8)

p <- profvis({fit.nm.rpt.loglinear(initial=initial, 
                                          known=known,
                                          data=train.df,
                                          weather=weather,
                                          include=include,
                                          print=TRUE,
                                          iter=1,
                                          scale=scale)
})
print(p)
#11740 ms before optimising
#11350 ms after unnaming parameter vectors in SEIR function
#10560 ms after storing derivatives in fewer steps