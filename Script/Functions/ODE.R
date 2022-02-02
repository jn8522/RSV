#functions related to the ODE solver
#the extended models here are the new ones that use four variables
library(deSolve)
library(tidyverse)
theme_set(theme_bw())

#defining the ODE model
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
  
  beta = beta0*(1+beta1*cos(2*pi*t/52.14286+phi))
  
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
  
  #pulling out the parameters that are being fitted
  nu = parms[1]
  p = parms[2]
  beta = betas[t]
  
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

#function that solves base ode model for given parameters and returns new cases per week
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
         0)

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
         0)
  
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

#for a logistic beta function
ode.solve.logistic <- function(variable, #variable parameters
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
         0)
  
  #calculating betas
  #setting unincluded parameters equal to 0
  cs <- rep(0, 4)
  L <- variable[3]
  c.0 <- variable[4]
  index=0
  for(i in 1:4){
    if(include[i]==TRUE){
      cs[i] = variable[i+4-index]
    }
    else{
      index = index+1
    }
  }
  betas = L/(1+exp(-c.0
                   -cs[1]*weather$precipitation
                   -cs[2]*weather$temperature
                   -cs[3]*weather$vapour
                   -cs[4]*weather$school))
  
  #making sure variables are in the right order
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

#function that computes the squared errors of a set of parameters theta wrt some data
#for the base model:
ssq.parameters <- function(theta,
                           data, #observed df
                           known){
  #unpacking the data
  fitted.df <- ode.solve(variable=theta, 
                         known=known, 
                         weeks=nrow(data),
                         noise=FALSE)
  #calculating error
  out <- mean((data$cases-fitted.df$cases)^2)
  return(out)
}

#function that computes the squared errors of a set of parameters theta wrt some data using the extended model
#for the loglinear model:
ssq.parameters.loglinear <- function(theta,
                           data, #observed df
                           known,
                           weather, #weather data covering same weeks as data
                           include){
  #calling the ODE solver with these parameters
  fitted.df <- ode.solve.loglinear(variable=theta, 
                         known=known,
                         noise=FALSE,
                         weather=weather,
                         include=include)
  #calculating error
  out <- mean((data$cases-fitted.df$cases)^2)
  return(out)
}
#for the logistic model:
ssq.parameters.logistic <- function(theta,
                                     data, #observed df
                                     known,
                                     weather, #weather data covering same weeks as data
                                     include){
  #calling the ODE solver with these parameters
  fitted.df <- ode.solve.logistic(variable=theta, 
                                   known=known,
                                   noise=FALSE,
                                   weather=weather,
                                   include=include)
  #calculating error
  out <- mean((data$cases-fitted.df$cases)^2)
  return(out)
}

#similar function that computes the (negative) log likelihood
#for the base model:
neglk.parameters <- function(theta,
                             data, #observed df
                             known){
  #solving ODE for the given parameters
  fitted.df <- ode.solve(variable=theta, 
                         known=known, 
                         weeks=nrow(data),
                         noise=FALSE) %>%
    filter(week>1)
  
  observed.df <- data %>% 
    filter(week>1)

  #calculating likelihood of observed data
  p = theta[5]
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
#for the extended model:
neglk.parameters.loglinear <- function(theta,
                             data, #observed df
                             known,
                             weather,
                             include){
  #solving ODE for the given parameters
  fitted.df <- ode.solve.loglinear(variable=theta, 
                         known=known,
                         noise=FALSE,
                         weather=weather,
                         include)%>%
    filter(week>1)
  
  observed.df <- data %>% filter(week>1)
  #calculating likelihood of observed data
  p = theta[2]
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

neglk.parameters.logistic <- function(theta,
                                       data, #observed df
                                       known,
                                       weather,
                                       include){
  #solving ODE for the given parameters
  fitted.df <- ode.solve.logistic(variable=theta, 
                                   known=known,
                                   noise=FALSE,
                                   weather=weather,
                                   include)%>%
    filter(week>1)
  
  observed.df <- data %>% filter(week>1)
  #calculating likelihood of observed data
  p = theta[2]
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

#overlay observed data with model output for parameter values theta
#for the base model:
plot.fit <- function(observed.df, #dataframe of observed data
                     variable, #variable parameters
                     known #fixed parameters
                     ) {
  fitted.df <- ode.solve(variable=variable,
                         known=known,
                         weeks=nrow(observed.df),
                         noise=FALSE)
  combined.df <- left_join(observed.df, fitted.df, 
                           by="week",
                           suffix=c(".observed", ".fitted"))%>%
    pivot_longer(!week, names_to = "source", values_to="value")%>%
      mutate(year = week/52.14286)
  end_year = max(combined.df$year) %>% ceiling()
  plot <- combined.df %>%
    ggplot(aes(x=year, y=value, color=source))+geom_line()+
    scale_x_continuous(breaks=seq(0,end_year,by=1))+
    scale_color_manual(values=c("red", "black"))+
    theme(legend.title=element_blank(), legend.position="bottom")+
    labs(y="Cases")
  return(plot)
}

#for the extended model
plot.fit.loglinear <- function(observed.df, #dataframe of observed data
                     variable, #variable parameters
                     known, #fixed parameters
                     weather,
                     include) {
  fitted.df <- ode.solve.loglinear(variable=variable,
                         known=known,
                         noise=FALSE,
                         weather=weather,
                         include=include)
  combined.df <- left_join(observed.df, fitted.df, 
                           by="week",
                           suffix=c(".observed", ".fitted")) %>%
    pivot_longer(!week, names_to = "source", values_to="value")%>%
    mutate(year = week/52.14286)
  end_year = max(combined.df$year) %>% ceiling()
  plot <- combined.df %>%
    ggplot(aes(x=year, y=value, color=source))+geom_line()+
    scale_x_continuous(breaks=seq(0,end_year,by=1))+
    scale_color_manual(values=c("red", "black"))+
    theme(legend.title=element_blank(), legend.position="bottom")+
    labs(y="Cases")
  return(plot)
}

plot.fit.logistic <- function(observed.df, #dataframe of observed data
                               variable, #variable parameters
                               known, #fixed parameters
                               weather,
                               include) {
  fitted.df <- ode.solve.logistic(variable=variable,
                             known=known,
                             noise=FALSE,
                             weather=weather,
                             include=include)
  combined.df <- left_join(observed.df, fitted.df, 
                           by="week",
                           suffix=c(".observed", ".fitted")) %>%
    pivot_longer(!week, names_to = "source", values_to="value")%>%
    mutate(year = week/52.14286)
  end_year = max(combined.df$year) %>% ceiling()
  plot <- combined.df %>%
    ggplot(aes(x=year, y=value, color=source))+geom_line()+
    scale_x_continuous(breaks=seq(0,end_year,by=1))+
    scale_color_manual(values=c("red", "black"))+
    theme(legend.title=element_blank(), legend.position="bottom")+
    labs(y="Cases")
  return(plot)
}

#function for fitting via repeated Nelder-Mead
#for the base model:
fit.nm.rpt <- function(initial, 
                       known, 
                       data, 
                       print, 
                       iter,
                       scale) {
  parm = initial
  n=0
 while(n<iter) {
    fit <- optim(par=parm,
                 fn = neglk.parameters,
                       control = list(trace=print,
                                      maxit=iter-n,
                                      parscale=scale),
                       known=known,
                       data=data)
    parm = fit$par
    n = n+fit$counts["function"]
    if(fit$convergence==0){
      break
    }
  }
  return(fit)
}

#for the extended model:
fit.nm.rpt.loglinear <- function(initial, 
                           known, 
                           data, 
                           weather, 
                           include, 
                           print, 
                           iter,
                           scale) {
  parm = initial
  n=0
  while(n<iter) {
    fit <- optim(par=parm,
                 fn = neglk.parameters.loglinear,
                 control = list(trace=print,
                                maxit=iter-n,
                                parscale=scale),
                 known=known,
                 data=data,
                 weather=weather,
                 include=include)
    parm = fit$par
    n = n+fit$counts["function"]
    if(fit$convergence==0){
      break
    }
  }
  return(fit)
}

fit.nm.rpt.logistic <- function(initial, 
                                 known, 
                                 data, 
                                 weather, 
                                 include, 
                                 print, 
                                 iter,
                                 scale) {
  parm = initial
  n=0
  while(n<iter) {
    fit <- optim(par=parm,
                 fn = neglk.parameters.logistic,
                 control = list(trace=print,
                                maxit=iter-n,
                                parscale=scale),
                 known=known,
                 data=data,
                 weather=weather,
                 include=include)
    parm = fit$par
    n = n+fit$counts["function"]
    if(fit$convergence==0){
      break
    }
  }
  return(fit)
}