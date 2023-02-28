source("src/fit.R")
source("src/ODE.R")
source("src/plot.R")
library(numDeriv)
library(dfoptim)
library(parallel)
#reading in the model fits
fit.base <- readRDS("models/fit.base.Rdata")
fit.weather <- readRDS("models/fit.weather.Rdata")

#reading in some other relevant data
RSV <- read.csv("data/processed/RSV.csv")
known <- c(read.csv("data/raw/known.values.csv")$value, RSV$cases[1])
weather.school <- read.csv("data/processed/weather.school.csv")%>%filter(week<=nrow(RSV))
include=c(TRUE, #c.p
          TRUE, #c.t
          FALSE, #c.v
          TRUE) #c.s

#WALD CIs
#Base model
hessian.base <- hessian(neglk.parameters.base,
                   fit.base$par,
                   data=RSV,
                   known=known,
                   i=0,
                   k=0)
write_rds(hessian.base,file="models/hessian.base.Rdata")
v.hat.base <- solve(hessian.base)
cis.wald.base <- data.frame(lower=fit.base$par-1.96*sqrt(diag(v.hat.base)),
                             upper=fit.base$par+1.96*sqrt(diag(v.hat.base)))
saveRDS(cis.wald.base,"data/statistics/cis.wald.base.RData")

#weather model
hessian.weather <- hessian(neglk.parameters.weather,
                           fit.weather$par,
                           data=RSV,
                           known=known,
                           weather=weather.school,
                           include=include,
                           i=0,
                           k=0)
write_rds(hessian.weather,file="models/hessian.weather.Rdata")
v.hat.weather <- solve(hessian.weather)
cis.wald.weather <- data.frame(lower=fit.weather$par-1.96*sqrt(diag(v.hat.weather)),
                            upper=fit.weather$par+1.96*sqrt(diag(v.hat.weather)))
saveRDS(cis.wald.weather,"data/statistics/cis.wald.weather.RData")

#PROFILE LIKELIHOODS
q = qchisq(0.95, 1)
l.star.base = neglk.parameters.base(fit.base$par, RSV, known,i=0,k=0)+q/2
l.star.weather = neglk.parameters.weather(fit.weather$par, RSV, known, weather.school, include,i=0,k=0)+q/2

#function that calculates log likelihood where parameters are split into beta (with associated index) and vector theta with other parameters.
#for the base model
neglk.wrapper.base <- function(par, beta, j, data, known){
  if(j==1){
    theta <- c(beta,par)
  }
  else if(j==11){
    theta <- c(par,beta)
  }
  else{
    theta <- c(par[1:(j-1)],beta,par[j:10])
  }
  return(min(neglk.parameters.base(theta, data, known,i=0,k=0),1e+20))
}

#for the weather model
neglk.wrapper.weather <- function(par, beta, j, data, known, weather, include){
  if(j==1){
    theta <- c(beta,par)
  }
  else if(j==12){
    theta <- c(par,beta)
  }
  else{
    theta <- c(par[1:(j-1)],beta,par[j:11])
  }
  return(min(neglk.parameters.weather(theta, data, known, weather,include,i=0,k=0),1e+20))
}

#function that calculates l-tilde(beta)
#i.e. if we're calculating the confidence interval for parameter beta,
#the function calculates the maximum likelihood of a particular value of beta by optimising over the other parameters
#the index argument j indicates which parameter beta is
#for the base model:
l.tilde.base <- function(theta.index, beta, j, data, known) {
  #bounds for optimiser (from which we need to omit beta)
  lower=c(0,0,0,0,0,0,0,0,0,0,0)
  upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1)
  
  if(j==1){
    par <- thetas.base[[theta.index]][2:11]
    lower.par <- lower[2:11]
    upper.par <- upper[2:11]
  }
  else if(j==11){
    par <- thetas.base[[theta.index]][1:10]
    lower.par <- lower[1:10]
    upper.par <- upper[1:10]
  }
  else{
    par <- c(thetas.base[[theta.index]][1:(j-1)], thetas.base[[theta.index]][(j+1):11])
    lower.par <- c(lower[1:(j-1)], lower[(j+1):11])
    upper.par <- c(upper[1:(j-1)], upper[(j+1):11])
  }
  fit <- hjkb(par,
              neglk.wrapper.base,
              lower=lower.par,
              upper=upper.par,
              control=list(info=TRUE),
              data=data,
              known=known,
              j=j,
              beta=beta)
  #checking that it definitely converged
  fit2 <- hjkb(fit$par,
              neglk.wrapper.base,
              lower=lower.par,
              upper=upper.par,
              control=list(info=TRUE),
              data=data,
              known=known,
              j=j,
              beta=beta)
  return(fit2)
}
#for the weather model:
l.tilde.weather <- function(theta.index, beta, j, data, known, weather, include) {
  #bounds for optimiser (from which we need to omit beta)
  lower=c(0,0,-Inf,-Inf,-Inf,-Inf,0,0,0,0,0,0)
  upper=c(Inf,1,Inf,Inf,Inf,Inf,1,1,1,1,1,1)
  
  if(j==1){
    par <- thetas.weather[[theta.index]][2:12]
    lower.par <- lower[2:12]
    upper.par <- upper[2:12]
  }
  else if(j==12){
    par <- thetas.weather[[theta.index]][1:10]
    lower.par <- lower[1:11]
    upper.par <- upper[1:11]
  }
  else{
    par <- c(thetas.weather[[theta.index]][1:(j-1)], thetas.weather[[theta.index]][(j+1):12])
    lower.par <- c(lower[1:(j-1)], lower[(j+1):12])
    upper.par <- c(upper[1:(j-1)], upper[(j+1):12])
  }
  fit <- hjkb(par,
              neglk.wrapper.weather,
              lower=lower.par,
              upper=upper.par,
              control=list(info=TRUE),
              data=data,
              known=known,
              weather=weather,
              include=include,
              j=j,
              beta=beta)
  #checking that it definitely converged
  fit2 <- hjkb(fit$par,
              neglk.wrapper.weather,
              lower=lower.par,
              upper=upper.par,
              control=list(info=TRUE),
              data=data,
              known=known,
              weather=weather,
              include=include,
              j=j,
              beta=beta)
  return(fit2)
}

#function that computes the error on a potential CI bound beta, i.e. how close beta is to q
error.beta.base <- function(beta,theta.index, j, data, known, l.star){
  l = l.tilde.base(theta.index, beta, j, data, known)
  if(j==1){
    thetas.base[[theta.index]] <<- c(beta,l$par)
  }
  else if(j==11){
    thetas.base[[theta.index]] <<- c(l$par,beta)
  }
  else{
    thetas.base[[theta.index]] <<- c(l$par[1:(j-1)],beta,l$par[j:10])
  }
  out = abs(l$value-l.star)
  print(out)
  print(thetas.base[[theta.index]])
  return(out)
}

error.beta.weather <- function(beta,theta.index, j, data, known, weather, include, l.star){
  l = l.tilde.weather(theta.index, beta, j, data, known, weather, include)
  if(j==1){
    thetas.weather[[theta.index]] <<- c(beta,l$par)
  }
  else if(j==12){
    thetas.weather[[theta.index]] <<- c(l$par,beta)
  }
  else{
    thetas.weather[[theta.index]] <<- c(l$par[1:(j-1)],beta,l$par[j:11])
  }
  out = abs(l$value-l.star)
  print(out)
  print(thetas.weather[[theta.index]])
  return(out)
}

#solving for beta_lower or beta_upper, as indicated by argument low = TRUE/FALSE
beta.solve.base <- function(j, low, beta.min, beta.max,theta.index, data, known, l.star) {
  if(low==TRUE){
    interval = c(beta.min,thetas.base[[theta.index]][j])
  }
  if(low==FALSE){
    interval=c(thetas.base[[theta.index]][j], beta.max)
  }
  fit <- optimise(error.beta.base,
                  interval=interval,
                  j=j,
                  data=data,
                  known=known,
                  l.star=l.star,
                  theta.index=theta.index)
  return(fit)
}

beta.solve.weather <- function(j, low, beta.min, beta.max,theta.index, data, known,weather,include, l.star) {
  if(low==TRUE){
    interval = c(beta.min,thetas.weather[[theta.index]][j])
  }
  if(low==FALSE){
    interval=c(thetas.weather[[theta.index]][j], beta.max)
  }
  fit <- optimise(error.beta.weather,
                  interval=interval,
                  j=j,
                  data=data,
                  known=known,
                  weather=weather,
                  include=include,
                  l.star=l.star,
                  theta.index=theta.index)
  return(fit)
}
#solving for the base model
j.base <- c(1,1,2,2,3,3,4,4,5,5)
theta.indexes.base <- 1:10
thetas.base <- list()
for (i in 1:10) {
  thetas.base[[i]] <- fit.base$par
}
lower.base <- rep(c(TRUE,FALSE),5)
beta.min.base <- c(0.8*fit.base$par["beta0"] %>% rep(2),
                   0.8*fit.base$par["beta1"] %>% rep(2),
                   0.8*fit.base$par["phi"] %>% rep(2),
                   0.5*fit.base$par["nu"] %>% rep(2),
                   0.5*fit.base$par["p"] %>% rep(2))
beta.max.base <- c(1.2*fit.base$par["beta0"] %>% rep(2),
                   1.2*fit.base$par["beta1"] %>% rep(2),
                   1.2*fit.base$par["phi"] %>% rep(2),
                   2*fit.base$par["nu"] %>% rep(2),
                   2*fit.base$par["p"] %>% rep(2))
cis.profile.base <- mcmapply(beta.solve.base,
                             mc.cores=10,
                           j=j.base,
                           theta.index=theta.indexes.base,
                           low=lower.base,
                           beta.min=beta.min.base,
                           beta.max=beta.max.base,
                           MoreArgs=list(data=RSV,
                                         known=known,
                                         l.star=l.star.base))
cis.profile.base["par"] <- list("beta0" %>% rep(2),
                              "beta1"%>% rep(2),
                              "phi" %>% rep(2),
                              "nu" %>% rep(2),
                              "p" %>% rep(2))
saveRDS(cis.profile.base,"data/statistics/cis.profile.base.RData")
#solving for the weather model
j.weather <- c(1,1,2,2,3,3,4,4,5,5,6,6)
theta.indexes <- 1:12
lower.weather <- rep(c(TRUE,FALSE),6)
beta.min.weather <- c(0.5*fit.weather$par["nu"] %>% rep(2),
                      0.5*fit.weather$par["p"] %>% rep(2),
                      0.5*fit.weather$par["c.0"] %>% rep(2),
                      0.5*fit.weather$par["c.p"] %>% rep(2),
                      0.5*fit.weather$par["c.t"] %>% rep(2),
                      0.5*fit.weather$par["c.s"] %>% rep(2))
beta.max.weather <- c(1.5*fit.weather$par["nu"] %>% rep(2),
                      1.5*fit.weather$par["p"] %>% rep(2),
                      1.5*fit.weather$par["c.0"] %>% rep(2),
                      1.5*fit.weather$par["c.p"] %>% rep(2),
                      1.5*fit.weather$par["c.t"] %>% rep(2),
                      1.5*fit.weather$par["c.s"] %>% rep(2))
thetas.weather <- list()
for (i in 1:12) {
  thetas.weather[[i]] <- fit.weather$par
}
cis.profile.weather <- mcmapply(beta.solve.weather,
                             mc.cores=12,
                             j=j.weather,
                             theta.index=theta.indexes,
                             low=lower.weather,
                             beta.min=beta.min.weather,
                             beta.max=beta.max.weather,
                             MoreArgs=list(data=RSV,
                                           known=known,
                                           l.star=l.star.weather,
                                           weather=weather.school,
                                           include=include))%>%
  data.frame()
colnames(cis.profile.weather) <- c("nu1",
                                   "nu2",
                                   "p1",
                                   "p2",
                                   "c.01",
                                   "c.02",
                                   "c.p1",
                                   "c.p2",
                                   "c.t1",
                                   "c.t2",
                                   "c.s1",
                                   "c.s2")

saveRDS(cis.profile.weather,"data/statistics/cis.profile.weather.RData")
