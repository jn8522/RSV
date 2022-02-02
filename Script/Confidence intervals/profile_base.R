#trying to find profile likelihoods for base model
library(numDeriv)
library(dfoptim)
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#reading in data
training_data()

#mle for base model
par.base <- read.csv("~/rsv-modelling/Output/Data/base.par.csv", header=FALSE)$V1

q = qchisq(0.95, 1)
l.star = neglk.parameters(par.base, train.df, known)+q/2

#function that calculates log likelihood where parameters are split into beta (with associated index) and vector theta with other parameters.
neglk.parameters.wrapper <- function(theta, beta, j, data, known){
  if(j==1){
    par <- c(beta,theta)
  }
  else if(j==11){
    par <- c(theta,beta)
  }
  else{
    par <- c(theta[1:(j-1)],beta,theta[j:10])
  }
  return(min(neglk.parameters(par, data, known),1e+20))
}

#function that calculates l-tilde(beta)
#i.e. if we're calculating the confidence interval for parameter beta,
#the function calculates the maximum likelihood of a particular value of beta by optimising over the other parameters
#the index argument j indicates which parameter beta is
l.tilde <- function(theta, beta, j, data, known) {
  #bounds for optimiser (from which we need to omit beta)
  lower=c(0,0,0,0,0,0,0,0,0,0,0)
  upper=c(Inf, 1, 2*pi, Inf, 1, 1,1,1,1,1,1)

  if(j==1){
    par <- theta[2:11]
    lower.par <- lower[2:11]
    upper.par <- upper[2:11]
  }
  else if(j==11){
    par <- theta[1:10]
    lower.par <- lower[1:10]
    upper.par <- upper[1:10]
  }
  else{
    par <- c(theta[1:(j-1)], theta[(j+1):11])
    lower.par <- c(lower[1:(j-1)], lower[(j+1):11])
    upper.par <- c(upper[1:(j-1)], upper[(j+1):11])
  }
  fit <- hjkb(par,
              neglk.parameters.wrapper,
              lower=lower.par,
              upper=upper.par,
              control=list(info=TRUE),
              j=j,
              beta=beta,
              data=data,
              known=known)
  return(fit)
}

#function that computes the error on a potential CI bound beta, i.e. how close beta is to q
error.beta <- function(beta, j, data, known){
  l = l.tilde(theta, beta, j, data, known)
  if(j==1){
    theta <<- c(beta,l$par)
  }
  else if(j==11){
    theta <<- c(l$par,beta)
  }
  else{
    theta <<- c(l$par[1:(j-1)],beta,l$par[j:10])
  }
  out = abs(l$value-l.star)
  return(out)
}

#solving for beta_lower or beta_upper, as indicated by argument low = TRUE/FALSE
beta_solve <- function(j, low, beta.min, beta.max) {
  if(low==TRUE){
    interval = c(beta.min,par.base[j])
  }
  if(low==FALSE){
    interval=c(par.base[j], beta.max)
  }
  fit <- optimise(error.beta,
              interval=interval,
              j=j,
              data=train.df,
              known=known)
}
theta=par.base
j=10
which="upper"
beta <- beta_solve(j=j,low=ifelse(which=="lower", TRUE, FALSE), beta.min=0, beta.max=0.3663191)
capture.output(beta, file=paste0("~/rsv-modelling/Output/Data/Confidence intervals/base_",j,"_",which,".txt"))