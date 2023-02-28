#functions that divide up a time series for the purpose of k-fold cross-validation
#these divide the training data into k subsamples and return just the ith (k.fold.test) or all but the ith (k.fold.train)
#setting k=0 makes training data=test data=whole time series
k.fold.train <- function(data,
                         i, #which subsample is currently serving as test data
                         k){ #total number of subsamples
  if(k==0){
    return(data)
  }else{
  index1 <- round((i-1)*1/k*nrow(data))+1
  index2 <-round(i*1/k*nrow(data))
  return(data %>% slice(-(index1:index2)))}
}
k.fold.test <- function(data,
                        i, #which subsample is currently serving as test data
                        k){ #total number of subsamples
  if(k==0){
    return(data)
  }else{
  index1 <- round((i-1)*1/k*nrow(data))+1
  index2 <-round(i*1/k*nrow(data))
  return(data %>% slice(index1:index2))}
}

#function that computes the mean squared errors of a set of parameters theta wrt some data
#the data is partitioned into k folds with the ith being reserved for test data and the remainder serving as training data
#the mse is calculated only over the test data
#for the base model:
mse.parameters.base <- function(theta,
                                data, #time series
                                known, #known parameters for the model
                                i,
                                k){
  #unpacking the data
  fitted.df <- ode.solve.base(variable=theta, 
                         known=known, 
                         weeks=nrow(data),
                         noise=FALSE)
  #calculating error
  out <- mean((k.fold.test(data,i,k)$cases-k.fold.test(fitted.df,i,k)$cases)^2)
  return(out)
}

#for the weather model:
mse.parameters.weather <- function(theta,
                                   data, #observed df
                                   known,
                                   weather, #weather data covering same weeks as data
                                   include,
                                   i,
                                   k){
  #calling the ODE solver with these parameters
  fitted.df <- ode.solve.weather(variable=theta, 
                                   known=known,
                                   noise=FALSE,
                                   weather=weather,
                                   include=include)
  #calculating error
  out <- mean((k.fold.test(data,i,k)$cases-k.fold.test(fitted.df,i,k)$cases)^2)
  return(out)
}

#similar function that computes the (negative) log likelihood
#the data is partitioned into k folds with the ith being reserved for test data and the remainder serving as training data
#the log likelihood is calculated only over the training data
#for the base model:
neglk.parameters.base <- function(theta,
                                  data, #observed df
                                  known,
                                  i,
                                  k){
  #solving ODE for the given parameters
  fitted.df <- ode.solve.base(variable=theta, 
                         known=known, 
                         weeks=nrow(data),
                         noise=FALSE) %>%
    filter(week>1)%>%
    k.fold.train(i,k)
  
  observed.df <- data %>% 
    filter(week>1)%>%
    k.fold.train(i,k)
  
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
#for the weather model:
neglk.parameters.weather <- function(theta,
                                     data, #observed df
                                     known,
                                     weather,
                                     include,
                                     i,
                                     k){
  #solving ODE for the given parameters
  fitted.df <- ode.solve.weather(variable=theta, 
                                   known=known,
                                   noise=FALSE,
                                   weather=weather,
                                   include)%>%
    filter(week>1)%>%
    k.fold.train(i,k)
  
  observed.df <- data%>%
    filter(week>1)%>%
    k.fold.train(i,k)
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