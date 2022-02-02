#calculating each model's forecast from the beginning of the study period through to 2021
#including a speculative reduction in transmission during Perth's 2020 lockdown
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
library(lubridate)

#usual preamble: training data, weather data, etc.
training_data()

#reading in test data
test.df <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/test.csv")
weather.train <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.train.csv")
weather.test <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.test.csv")
weather.forecast <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.forecast.csv")
weather <- union(weather.train, weather.test) %>%
  union(weather.forecast) %>%
  mutate(precipitation = precipitation/mean(weather.train$precipitation),
         evaporation = evaporation/mean(weather.train$evaporation),
         temperature = temperature/mean(weather.train$temperature),
         dewpoint = dewpoint/mean(weather.train$dewpoint),
         humidity = humidity/mean(weather.train$humidity))

#reading in the fitted parameters
par.base = read.csv("~/rsv-modelling/Output/Data/base.par.csv", header=FALSE)$V1
par.loglinear <- read.csv("~/rsv-modelling/Output/Data/Fits/Backwards selection/Loglinear round 1/Fit_loglinear_v_par.csv")$x
include=c(TRUE,TRUE,FALSE,TRUE)

#Stage 1 restrictions: 21 March
#Stage 2 restrictions: 23 March
interval(dmy("1.1.2000"),dmy("23.3.2020"))%/% weeks(1)+1
#Stage 3 restrictions: 29 March
#Phase 1 easing: 26 April
interval(dmy("1.1.2000"),dmy("26.4.2020"))%/% weeks(1)+1
#Schools resume: 30 April
interval(dmy("1.1.2000"),dmy("30.4.2020"))%/% weeks(1)+1
#Phase 2 easing: 16 May
interval(dmy("1.1.2000"),dmy("16.5.2020"))%/% weeks(1)+1
#Phase 3 easing: 5 June
interval(dmy("1.1.2000"),dmy("5.6.2020"))%/% weeks(1)+1

#redefining all the ode functions, adding a lockdown of variable efficacy from 21 March to 5 June
seir.model <- function(t, y, parms, known) {
  #pulling out the state
  S1 = y[1]
  E1 = y[2]
  I1 = y[3]
  R1 = y[4]
  S2 = y[5]
  E2 = y[6]
  I2 = y[7]
  R2 = y[8]
  
  #pulling out the parameters that are being fitted
  beta0 = parms[1]
  beta1 = parms[2]
  phi = parms[3]
  nu = parms[4]
  p = parms[5]
  
  #pulling out the known parameters
  mu = known[1]
  eta1 = known[2]
  eta2 = known[3]
  alpha = known[4]
  sigma = known[5]
  gamma = known[6]
  delta = known[7]
  n1 = known[8]
  n2 = known[9]
  efficacy = known[11]
  
  beta = beta0*(1+beta1*cos(2*pi*t/52.14286+phi))
  if(1056<=t&t<=1066){
    beta=beta*(1-efficacy/100)
  }
  
  #model equations
  dS1dt = mu*n2-beta*S1*((I1+alpha*I2)/(n1+n2))+nu*R1-eta1*S1
  dE1dt = beta*S1*((I1+alpha*I2)/(n1+n2))-sigma*E1-eta1*E1
  dI1dt = sigma*E1-gamma*I1-eta1*I1
  dR1dt = gamma*I1-nu*R1-eta1*R1
  dS2dt = eta1*S1-delta*beta*S2*((I1+alpha*I2)/(n1+n2))+nu*R2-eta2*S2
  dE2dt = eta1*E1+delta*beta*S2*((I1+alpha*I2)/(n1+n2))-sigma*E2-eta2*E2
  dI2dt = eta1*I1+sigma*E2-gamma*I2-eta2*I2
  dR2dt = eta1*R1+gamma*I2-nu*R2-eta2*R2
  
  #new cases in the 0-2 age group
  new = sigma*E1
  
  dydt = c(dS1dt, dE1dt, dI1dt, dR1dt, dS2dt, dE2dt, dI2dt, dR2dt, new)
  
  #return
  list(dydt)
}

#defining the extended ODE model
#weather is a data frame
seir.model.ext <- function(t, y, parms, known, betas) {
  #pulling out the state
  S1 = y[1]
  E1 = y[2]
  I1 = y[3]
  R1 = y[4]
  S2 = y[5]
  E2 = y[6]
  I2 = y[7]
  R2 = y[8]
  
  #pulling out the known parameters
  mu = known[1]
  eta1 = known[2]
  eta2 = known[3]
  alpha = known[4]
  sigma = known[5]
  gamma = known[6]
  delta = known[7]
  n1 = known[8]
  n2 = known[9]
  efficacy = known[11]
  
  #pulling out the parameters that are being fitted
  nu = parms[1]
  p = parms[2]
  beta = betas[t]
  
  if(1056<=t&t<=1066){
    beta=beta*(1-efficacy/100)
  }
  
  list(c(mu*n2-beta*S1*((I1+alpha*I2)/(n1+n2))+nu*R1-eta1*S1, #dS1/dt
         beta*S1*((I1+alpha*I2)/(n1+n2))-sigma*E1-eta1*E1, #dE1/dt
         sigma*E1-gamma*I1-eta1*I1, #dI1/dt
         gamma*I1-nu*R1-eta1*R1, #dR1/dt
         eta1*S1-delta*beta*S2*((I1+alpha*I2)/(n1+n2))+nu*R2-eta2*S2, #dS2/dt
         eta1*E1+delta*beta*S2*((I1+alpha*I2)/(n1+n2))-sigma*E2-eta2*E2, #dE2/dt
         eta1*I1+sigma*E2-gamma*I2-eta2*I2, #dI2/dt
         eta1*R1+gamma*I2-nu*R2-eta2*R2, #dR2/dt
         sigma*E1)) #new cases
}

ode.solve <- function(variable, #variable parameters
                      known, #fixed parameters
                      weeks, #how many weeks to solve for
                      noise #true/false, whether to add random noise
) {
  
  #determining initial conditions based on provided ratios
  n1 = known[8]
  n2 = known[9]
  ra1 = variable[6]
  ra2 = variable[7]
  ra3 = variable[8]
  ra4 = variable[9]
  ra5 = variable[10]
  ra6 = variable[11]
  s1r1 = (1-ra1)*n1
  e1r1 = ra1*n1
  s1 = (1-ra2)*s1r1
  r1 = ra2*s1r1
  e1 = (1-ra3)*e1r1
  i1 = ra3*e1r1
  s2r2 = (1-ra4)*n2
  e2r2 = ra4*n2
  s2 = (1-ra5)*s2r2
  r2 = ra5*s2r2
  e2 = (1-ra6)*e2r2
  i2 = ra6*e2r2
  p=variable[5]
  y0 = c(s1,e1,i1,r1,s2,e2,i2,r2,
         known[10]/p)
  
  #solving ODE and doing some data wrangling
  ode.output <- ode(func=seir.model, 
                    y=y0, 
                    times=1:weeks,
                    parms=variable,
                    known=known) %>%
    as.data.frame()
  
  names(ode.output) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new")
  
  #working out weekly new cases from cumulative new cases
  cases <- c(ode.output$new[1],(ode.output$new[2:weeks]-ode.output$new[1:weeks-1]))
  #adding random noise (if applicable)
  
  if(noise==TRUE){
    n=length(cases)
    noise.means <- cases*p
    noise.stdevs <- sqrt(cases*p*(1-p))
    cases.noise <- rnorm(n, noise.means, noise.stdevs)
    out <- data.frame(week = 1:length(cases),
                      cases = cases.noise)
  }
  else{
    out <- data.frame(week = 1:length(cases),
                      cases = cases*p)
  }
  return(out)
}

#functions that solve extended ode model for given parameters and returns new cases per week
#solves for the weeks covered by the weather data (no need to specify number of weeks)
#for a loglinear beta:
ode.solve.loglinear <- function(variable, #variable parameters
                                known, #fixed parameters
                                noise, #true/false, whether to add random noise
                                weather, #weather data
                                include) #named vector saying which variables to include
{
  #times to solve ODE for
  weeks=nrow(weather)
  
  #determining initial conditions based on provided ratios
  n1 = known[8]
  n2 = known[9]
  ra1 = variable[length(variable)-5]
  ra2 = variable[length(variable)-4]
  ra3 = variable[length(variable)-3]
  ra4 = variable[length(variable)-2]
  ra5 = variable[length(variable)-1]
  ra6 = variable[length(variable)]
  s1r1 = (1-ra1)*n1
  e1r1 = ra1*n1
  s1 = (1-ra2)*s1r1
  r1 = ra2*s1r1
  e1 = (1-ra3)*e1r1
  i1 = ra3*e1r1
  s2r2 = (1-ra4)*n2
  e2r2 = ra4*n2
  s2 = (1-ra5)*s2r2
  r2 = ra5*s2r2
  e2 = (1-ra6)*e2r2
  i2 = ra6*e2r2
  p = variable[2]
  y0 = c(s1,e1,i1,r1,s2,e2,i2,r2,
         known[10]/p)
  
  #calculating betas
  #setting unincluded parameters equal to 0
  cs <- rep(0, 4)
  c.0 <- variable[3]
  index=0
  for(i in 1:4){
    if(include[i]==TRUE){
      cs[i] = variable[i+3-index]
    }
    else{
      index = index+1
    }
  }
  betas = exp(c.0
              +cs[1]*weather$precipitation
              +cs[2]*weather$temperature
              +cs[3]*weather$vapour
              +cs[4]*weather$school)
  
  #solving ODE and doing some data wrangling
  ode.output <- ode(func=seir.model.ext, 
                    y=y0, 
                    times=1:weeks,
                    parms=variable,
                    known=known,
                    betas=betas) %>%
    as.data.frame()
  
  names(ode.output) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new")
  
  
  #working out weekly new cases from cumulative new cases
  cases <- c(ode.output$new[1],(ode.output$new[2:weeks]-ode.output$new[1:weeks-1]))
  
  #adding random noise (if applicable)
  if(noise==TRUE){
    n=length(cases)
    noise.means <- cases*p
    noise.stdevs <- sqrt(cases*p*(1-p))
    cases.noise <- rnorm(n, noise.means, noise.stdevs)
    out <- data.frame(week = 1:length(cases),
                      cases = cases.noise)
  }
  else{
    out <- data.frame(week = 1:length(cases),
                      cases = cases*p)
  }
  return(out)
}

#Forecasts
forecast.loglinear <- ode.solve.loglinear(variable=par.loglinear,
                                            known=c(known, 75),
                                            noise=FALSE,
                                            weather=weather,
                                            include=include)

forecast.base <- ode.solve(variable=par.base,
                             weeks=nrow(weather),
                             known=c(known, 75),
                             noise=FALSE)

library(RColorBrewer)
study.endyear = weather.test$week[nrow(weather.test)]/(365.25/7)+2000
loglinear.plot <- forecast.loglinear %>% 
  mutate(year=week/(365.25/7)+2000) %>%
  ggplot(aes(x=year, y=cases))+geom_line()+
  scale_x_continuous(breaks=seq(2000,2020,by=5))+
  labs(x="", y="Cases")+
  ggtitle("Extended model")+
  scale_color_brewer(palette="Dark2", name="Lockdown efficacy")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

base.plot <- forecast.base %>% 
  mutate(year=week/(365.25/7)+2000) %>%
  ggplot(aes(x=year, y=cases))+geom_line()+
  scale_x_continuous(breaks=seq(2000,2020,by=5))+
  labs(x="", y="Cases")+
  ggtitle("Sinusoidal model")+
  scale_color_brewer(palette="Dark2", name="Lockdown efficacy")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

library(ggpubr)
arranged.plot <- ggarrange(base.plot, loglinear.plot, nrow=2, ncol=1, common.legend=TRUE, legend="bottom")
ggsave(filename="forecast.full.tif",
         device="tiff",
         plot=arranged.plot,
         path="~/rsv-modelling/Output/Plots",
         width=15,
         height=15,
         units="cm",
         dpi=300)

