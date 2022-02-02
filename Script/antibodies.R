#calculating proportion who've ever been infected

#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#combining training data for training and test period
weather.train <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.train.csv")
weather.test <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.test.csv")
weather <- union(weather.train, weather.test) %>%
  mutate(precipitation = precipitation/mean(weather.train$precipitation),
         evaporation = evaporation/mean(weather.train$evaporation),
         temperature = temperature/mean(weather.train$temperature),
         dewpoint = dewpoint/mean(weather.train$dewpoint),
         humidity = humidity/mean(weather.train$humidity))

#reading in parameters
par.base <- read.csv("~/rsv-modelling/Output/Data/base.par.csv", header=FALSE)$V1
par.loglinear <- read.csv("~/rsv-modelling/Output/Data/Fits/Backwards selection/Loglinear round 1/Fit_loglinear_v_par.csv")$x
include=c(TRUE,TRUE,FALSE,TRUE)

#from ode.R
#introducing new variable A1 for number of 0--2s with antibodies
#for the base model:
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
  SA = y[10]
  EA = y[11]
  IA = y[12]
  
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
  dS1dt = mu*n2-beta*S1*(((I1+IA)+alpha*I2)/(n1+n2))-eta1*S1
  dE1dt = beta*S1*(((I1+IA)+alpha*I2)/(n1+n2))-sigma*E1-eta1*E1
  dI1dt = sigma*E1-gamma*I1-eta1*I1
  dR1dt = gamma*I1-nu*R1-eta1*R1+gamma*IA
  dSAdt = nu*R1-beta*SA*(((I1+IA)+alpha*I2)/(n1+n2))-eta1*SA
  dEAdt = beta*SA*(((I1+IA)+alpha*I2)/(n1+n2))-sigma*EA-eta1*EA
  dIAdt = sigma*EA-gamma*IA-eta1*IA
  dS2dt = eta1*S1+eta1*SA-delta*beta*S2*(((I1+IA)+alpha*I2)/(n1+n2))+nu*R2-eta2*S2
  dE2dt = eta1*E1+eta1*EA+delta*beta*S2*(((I1+IA)+alpha*I2)/(n1+n2))-sigma*E2-eta2*E2
  dI2dt = eta1*I1+eta1*IA+sigma*E2-gamma*I2-eta2*I2
  dR2dt = eta1*R1+gamma*I2-nu*R2-eta2*R2

  
  #new cases in the 0-2 age group
  new = sigma*E1
  
  dydt = c(dS1dt, dE1dt, dI1dt, dR1dt, dS2dt, dE2dt, dI2dt, dR2dt, new, dSAdt, dEAdt, dIAdt)
  
  #return
  list(dydt)
}

#solving
n1 = known[8]
n2 = known[9]
ra1 = par.base[6]
ra2 = par.base[7]
ra3 = par.base[8]
ra4 = par.base[9]
ra5 = par.base[10]
ra6 = par.base[11]
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
p=par.base[5]
y0 = c(s1,e1,i1,r1,s2,e2,i2,r2,
       known[10]/p, 
       0,0,0)

ode.output.base <- ode(func=seir.model, 
                       y=y0, 
                       times=1:nrow(weather),
                       parms=par.base,
                       known=known) %>%
  as.data.frame()
names(ode.output.base) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new", "SA","EA","IA")
ode.output.base <- ode.output.base %>%
  mutate(antibodies=SA+EA+IA+R1,
         n1=S1+E1+I1+R1+SA+EA+IA)

#for the extended model:
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
  SA = y[10]
  EA = y[11]
  IA = y[12]
  
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
  
  dS1dt = mu*n2-beta*S1*(((I1+IA)+alpha*I2)/(n1+n2))-eta1*S1
  dE1dt = beta*S1*(((I1+IA)+alpha*I2)/(n1+n2))-sigma*E1-eta1*E1
  dI1dt = sigma*E1-gamma*I1-eta1*I1
  dR1dt = gamma*I1-nu*R1-eta1*R1+gamma*IA
  dSAdt = nu*R1-beta*SA*(((I1+IA)+alpha*I2)/(n1+n2))-eta1*SA
  dEAdt = beta*SA*(((I1+IA)+alpha*I2)/(n1+n2))-sigma*EA-eta1*EA
  dIAdt = sigma*EA-gamma*IA-eta1*IA
  dS2dt = eta1*S1+eta1*SA-delta*beta*S2*(((I1+IA)+alpha*I2)/(n1+n2))+nu*R2-eta2*S2
  dE2dt = eta1*E1+eta1*EA+delta*beta*S2*(((I1+IA)+alpha*I2)/(n1+n2))-sigma*E2-eta2*E2
  dI2dt = eta1*I1+eta1*IA+sigma*E2-gamma*I2-eta2*I2
  dR2dt = eta1*R1+gamma*I2-nu*R2-eta2*R2
  
  
  #new cases in the 0-2 age group
  new = sigma*E1
  
  dydt = c(dS1dt, dE1dt, dI1dt, dR1dt, dS2dt, dE2dt, dI2dt, dR2dt, new, dSAdt, dEAdt, dIAdt)
  
  #return
  list(dydt)}

#solving
n1 = known[8]
n2 = known[9]
ra1 = par.loglinear[length(par.loglinear)-5]
ra2 = par.loglinear[length(par.loglinear)-4]
ra3 = par.loglinear[length(par.loglinear)-3]
ra4 = par.loglinear[length(par.loglinear)-2]
ra5 = par.loglinear[length(par.loglinear)-1]
ra6 = par.loglinear[length(par.loglinear)]
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
p = par.loglinear[2]
y0 = c(s1,e1,i1,r1,s2,e2,i2,r2,
       known[10]/p, 
       0,0,0)
cs <- rep(0, 4)
c.0 <- par.loglinear[3]
index=0
for(i in 1:4){
  if(include[i]==TRUE){
    cs[i] = par.loglinear[i+3-index]
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

ode.output.ext <- ode(func=seir.model.ext, 
                      y=y0, 
                      times=1:nrow(weather),
                      parms=par.loglinear,
                      known=known,
                      betas=betas) %>%
  as.data.frame()

names(ode.output.ext) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new", "SA","EA","IA")
ode.output.ext <- ode.output.ext %>%
  mutate(antibodies=SA+EA+IA+R1,
         n1=S1+E1+I1+R1+SA+EA+IA)

antibodies.df <- data.frame(year = ode.output.base$week/(365.25/7),
                         base = ode.output.base$antibodies/ode.output.base$n1,
                         ext = ode.output.ext$antibodies/ode.output.ext$n1)
antibodies.plot <- antibodies.df %>%
  pivot_longer(!year, names_to="Model", values_to="value")%>%
  ggplot(aes(x=year, y=value, colour=Model))+geom_line()
ggsave("Output/Plots/antibodies.plot.png",antibodies.plot)

mean(filter(antibodies.df, year>=2.5)$base)
mean(filter(antibodies.df, year>=2.5)$ext)

dist <- data.frame(draw=rexp(10000, rate=1/2))%>% 
  ggplot(aes(x=draw))+
  geom_density()+
  labs(x="Age (years)",y="")
ggsave("Output/Plots/dist.plot.png",dist)

mean(dist$draw)
