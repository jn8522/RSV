#this script calculates summary statistics assessing the model fits
source("src/fit.R")
source("src/ODE.R")
source("src/plot.R")

#reading in the model fits
fits.base <- readRDS("models/fits.base.CV.Rdata")
fits.weather <- readRDS("models/fits.weather.CV.Rdata")

#reading in some other relevant data
RSV <- read.csv("data/processed/RSV.csv")
known <- c(read.csv("data/raw/known.values.csv")$value, RSV$cases[1])
weather.school <- read.csv("data/processed/weather.school.csv")%>%filter(week<=nrow(RSV))
include=c(TRUE, #c.p
          TRUE, #c.t
          FALSE, #c.v
          TRUE) #c.s

#calculating the MSEs
MSE <- data.frame(i=1:5,
                  base=c(mse.parameters.base(fits.base[1,1]$par, data=RSV, known=known, i=1, k=5),
                         mse.parameters.base(fits.base[1,2]$par, data=RSV, known=known, i=2, k=5),
                         mse.parameters.base(fits.base[1,3]$par, data=RSV, known=known, i=3, k=5),
                         mse.parameters.base(fits.base[1,4]$par, data=RSV, known=known, i=4, k=5),
                         mse.parameters.base(fits.base[1,5]$par, data=RSV, known=known, i=5, k=5)),
                  weather=c(mse.parameters.weather(fits.weather[1,1]$par, data=RSV, known=known, weather=weather.school, include=include, i=1, k=5),
                            mse.parameters.weather(fits.weather[1,2]$par, data=RSV, known=known, weather=weather.school, include=include, i=2, k=5),
                            mse.parameters.weather(fits.weather[1,3]$par, data=RSV, known=known, weather=weather.school, include=include, i=3, k=5),
                            mse.parameters.weather(fits.weather[1,4]$par, data=RSV, known=known, weather=weather.school, include=include, i=4, k=5),
                            mse.parameters.weather(fits.weather[1,5]$par, data=RSV, known=known, weather=weather.school, include=include, i=5, k=5)))
avg.MSE <- c(base=mean(MSE$base), weather=mean(MSE$weather))

#also checking the MSE for each model over the entire timeseries when fitted to the entire timeseries
fit.base <- readRDS("models/fit.base.Rdata")
fit.weather <- readRDS("models/fit.weather.Rdata")
mse.parameters.base(fit.base$par, data=RSV,known=known,i=0,k=0)
mse.parameters.weather(fit.weather$par, data=RSV,known=known,weather=weather.school,include=include,i=0,k=0)

#Hierarchical partitioning
library(hier.part)
fit.weather <- readRDS("models/fit.weather.Rdata")
predictors <- data.frame(precipitation=weather.school$precipitation,
                         temperature=weather.school$temperature,
                         school=weather.school$school)
beta=exp(fit.weather$par["c.0"]
         +fit.weather$par["c.p"]*weather.school$precipitation
         +fit.weather$par["c.t"]*weather.school$temperature
         +fit.weather$par["c.s"]*weather.school$school)
hier <- hier.part(beta,predictors,family="poisson",gof="RMSPE",barplot=FALSE)

#checking whether sensitivity analysis changes the variance importance rankings of our variables
fits.weather <- readRDS("models/fits.sensitivity.weather.Rdata")
hiers<-list()
for (i in 1:20) {
  beta=exp(fits.weather[1,i]$par["c.0"]
           +fits.weather[1,i]$par["c.p"]*weather.school$precipitation
           +fits.weather[1,i]$par["c.t"]*weather.school$temperature
           +fits.weather[1,i]$par["c.s"]*weather.school$school)
  hiers[[i]] <- hier.part(beta,predictors,family="poisson",gof="RMSPE",barplot=FALSE)
}
