source("src/fit.R")
source("src/plot.R")
source("src/ODE.R")
library(ggpubr)
library(lubridate)
theme_set(theme_bw())

#reading in the model fits
fit.base <- readRDS("models/fit.base.Rdata")
fit.weather <- readRDS("models/fit.weather.Rdata")

#reading in some other relevant data
RSV <- read.csv("data/processed/RSV.csv")
known <- c(read.csv("data/raw/known.values.csv")$value, RSV$cases[1])
weather.school <- read.csv("data/processed/weather.school.csv")
include=c(TRUE, #c.p
          TRUE, #c.t
          FALSE, #c.v
          TRUE) #c.s

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

#redefining the ode functions, adding a lockdown of variable efficacy from 23 March to 5 June
seir.model.base <- function(t, y, parms, known) {
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
  n1 = S1+E1+I1+R1
  n2 = S2+E2+I2+R2
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
seir.model.weather <- function(t, y, parms, known, betas) {
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
  efficacy = known[11]
  if(1056<=t&t<=1066){
    beta=beta*(1-efficacy/100)
  }
  
  #pulling out the known parameters
  mu = known[1]
  eta1 = known[2]
  eta2 = known[3]
  alpha = known[4]
  sigma = known[5]
  gamma = known[6]
  delta = known[7]
  n1=S1+E1+I1+R1
  n2=S2+E2+I2+R2
  
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
ode.solve.base <- function(variable, #variable parameters
                           known, #fixed parameters
                           weeks, #how many weeks to solve for
                           noise #true/false, whether to add random noise
) {
  
  #determining initial conditions based on provided ratios
  #these are:
  #ra1: (E1+I1)/(S1+E1+I1+R1)
  #ra2: R1/(S1+R1)
  #ra3: I1/(E1+I1)
  #ra4: (E2+I2)/(S2+E2+I2+R2)
  #ra5: R2/(S2+R2)
  #ra6: I2/(E2+I2)
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
  ode.output <- ode(func=seir.model.base, 
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
                      cases = cases*p,
                      S1=ode.output$S1,
                      E1=ode.output$E1,
                      I1=ode.output$I1,
                      R1=ode.output$R1,
                      S2=ode.output$S2,
                      E2=ode.output$E2,
                      I2=ode.output$I2,
                      R2=ode.output$R2)
  }
  return(out)
}

#functions that solve extended ode model for given parameters and returns new cases per week
#solves for the weeks covered by the weather data (no need to specify number of weeks)
ode.solve.weather <- function(variable, #variable parameters
                              known, #fixed parameters
                              noise, #true/false, whether to add random noise
                              weather, #weather data
                              include) #named vector saying which variables to include
{
  #times to solve ODE for
  weeks=nrow(weather.school)
  
  #determining initial conditions based on provided ratios
  #these are:
  #ra1: (E1+I1)/(S1+E1+I1+R1)
  #ra2: R1/(S1+R1)
  #ra3: I1/(E1+I1)
  #ra4: (E2+I2)/(S2+E2+I2+R2)
  #ra5: R2/(S2+R2)
  #ra6: I2/(E2+I2)
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
  ode.output <- ode(func=seir.model.weather, 
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
                      cases = cases*p,
                      S1=ode.output$S1,
                      E1=ode.output$E1,
                      I1=ode.output$I1,
                      R1=ode.output$R1,
                      S2=ode.output$S2,
                      E2=ode.output$E2,
                      I2=ode.output$I2,
                      R2=ode.output$R2)
  }
  return(out)
}

#producing the forecasts
x100.base <- ode.solve.base(fit.base$par,
                            known=c(known,100),
                            weeks=nrow(weather.school),
                            noise=FALSE)
x75.base <- ode.solve.base(fit.base$par,
                           known=c(known,75),
                           weeks=nrow(weather.school),
                           noise=FALSE)
x50.base <- ode.solve.base(fit.base$par,
                           known=c(known,50),
                           weeks=nrow(weather.school),
                           noise=FALSE)
forecasts.base <- data.frame(week=weather.school$week,
                             x100=x100.base$cases,
                             x75=x75.base$cases,
                             x50=x50.base$cases)%>%
  pivot_longer(cols=!week,names_to="Efficacy",values_to="cases")%>%
  mutate(Efficacy=factor(Efficacy, levels=c("x50","x75","x100"), labels=c("50%","75%","100%")))

x100.weather <- ode.solve.weather(fit.weather$par,
                                 known=c(known,100),
                                 noise=FALSE,
                                 weather=weather.school,
                                 include=include)
x75.weather <- ode.solve.weather(fit.weather$par,
                                 known=c(known,75),
                                 noise=FALSE,
                                 weather=weather.school,
                                 include=include)
x50.weather <- ode.solve.weather(fit.weather$par,
                                 known=c(known,50),
                                 noise=FALSE,
                                 weather=weather.school,
                                 include=include)
forecasts.weather <- data.frame(week=weather.school$week,
                                x100=x100.weather$cases,
                                x75=x75.weather$cases,
                                x50=x50.weather$cases)%>%
  pivot_longer(cols=!week,names_to="Efficacy",values_to="cases")%>%
  mutate(Efficacy=factor(Efficacy, levels=c("x50","x75","x100"), labels=c("50%","75%","100%")))

#figure4
library(RColorBrewer)
study.endyear = weather.school$week[nrow(weather.school)]/(365.25/7)+2000
weather.plot <- forecasts.weather %>% 
  mutate(year=week/(365.25/7)+2000) %>%
  filter(year<2020)%>%
  ggplot(aes(x=year, y=cases))+geom_line()+
  scale_x_continuous(breaks=seq(2000,2020,by=5))+
  labs(x="", y="Weekly hospitalizations")+
  ggtitle("Weather model")+
  scale_color_brewer(palette="Dark2", name="Lockdown efficacy")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

base.plot <- forecasts.base %>% 
  mutate(year=week/(365.25/7)+2000) %>%
  filter(year<2020)%>%
  ggplot(aes(x=year, y=cases))+geom_line()+
  scale_x_continuous(breaks=seq(2000,2020,by=5))+
  labs(x="", y="Weekly hospitalizations")+
  ggtitle("Sinusoidal model")+
  scale_color_brewer(palette="Dark2", name="Lockdown efficacy")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

arranged.plot <- ggarrange(base.plot, weather.plot, nrow=2, ncol=1, common.legend=TRUE, legend="bottom")
ggsave(filename="Figure4.tif",
       device="tiff",
       plot=arranged.plot,
       path="plots",
       width=15,
       height=15,
       units="cm",
       dpi=300)

#figure 5
weather.plot2 <- forecasts.weather %>% 
  mutate(year=as.Date(date_decimal(week/(365.25/7)+2000))) %>%
  filter(year>=dmy("1/1/2020"))%>%
  ggplot(aes(x=year, y=cases, color=Efficacy))+geom_line()+
  scale_x_date(date_breaks="3 months", date_labels="%b %Y", minor_breaks="1 month")+
  geom_vline(xintercept=dmy("23/3/2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("27/4/2020"), y=60, label="Lockdown", size=3)+
  geom_vline(xintercept=dmy("1.6.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.7.2020"), y=60, label="Winter", size=3)+
  geom_vline(xintercept=dmy("1.9.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.10.2020"), y=60, label="Spring", size=3)+
  geom_vline(xintercept=dmy("1.12.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.1.2021"), y=60, label="Summer", size=3)+
  geom_vline(xintercept=dmy("1.3.2021"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.4.2021"), y=60, label="Autumn", size=3)+
  geom_vline(xintercept=dmy("1.6.2021"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.7.2021"), y=60, label="Winter", size=3)+
  labs(x="", y="Weekly hospitalizations")+
  ggtitle("Weather model")+
  scale_color_brewer(palette="Dark2", name="Lockdown efficacy")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

base.plot2 <- forecasts.base %>% 
  mutate(year=as.Date(date_decimal(week/(365.25/7)+2000))) %>%
  filter(year>=dmy("1/1/2020"))%>%
  ggplot(aes(x=year, y=cases, color=Efficacy))+geom_line()+
  scale_x_date(date_breaks="3 months", date_labels="%b %Y", minor_breaks="1 month")+
  geom_vline(xintercept=dmy("23/3/2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("27/4/2020"), y=60, label="Lockdown", size=3)+
  geom_vline(xintercept=dmy("1.6.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.7.2020"), y=60, label="Winter", size=3)+
  geom_vline(xintercept=dmy("1.9.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.10.2020"), y=60, label="Spring", size=3)+
  geom_vline(xintercept=dmy("1.12.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.1.2021"), y=60, label="Summer", size=3)+
  geom_vline(xintercept=dmy("1.3.2021"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.4.2021"), y=60, label="Autumn", size=3)+
  geom_vline(xintercept=dmy("1.6.2021"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.7.2021"), y=60, label="Winter", size=3)+
  labs(x="", y="Weekly hospitalizations")+
  ggtitle("Sinusoidal model")+
  scale_color_brewer(palette="Dark2", name="Lockdown efficacy")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

arranged.plot2 <- ggarrange(base.plot2, weather.plot2, nrow=2, ncol=1, common.legend=TRUE, legend="bottom")
ggsave(filename="Figure5.tif",
       device="tiff",
       plot=arranged.plot2,
       path="plots",
       width=15,
       height=15,
       units="cm",
       dpi=300)

#supplementary figure
#plotting R_eff under the different lockdown scenarios
r0 <- function(data, beta, known){
  e1 <- known[2]
  e2 <- known[3]
  a <- known[4]
  s <- known[5]
  g <- known[6]
  d <- known[7]
  n1 <- data$S1+data$E1+data$I1+data$R1
  n2 <- data$S2+data$E2+data$I2+data$R2
  out <- (beta*s*(a*e1^2*d*n2 + a*e1^2*n1 + a*e1*d*g*n2 + a*e1*d*s*n2 + a*e1*g*n1 + a*e1*s*n1 + a*e1*e2*n1 + a*d*g*s*n2 + g*s*n1 + g*e2*n1 + s*e2*n1 + e2^2*n1))/((e1 + g)*(e1 + s)*(g + e2)*(s + e2)*(n2 + n1))
  return(out)
}

rt <- function(data, beta, known){
  s1 = data$S1
  s2 = data$S2
  n1=data$S1+data$E1+data$I1+data$R1
  n2=data$S2+data$E2+data$I2+data$R2
  e1 <- known[2]
  e2 <- known[3]
  a <- known[4]
  s <- known[5]
  g <- known[6]
  d <- known[7]
  n1 <- known[8]
  n2 <- known[9]
  out = (beta*s*(a*e1^2*d*s2 + a*e1^2*s1 + a*e1*d*g*s2 + a*e1*d*s*s2 + a*e1*g*s1 + a*e1*s1*s + a*e1*s1*e2 + a*d*g*s*s2 + g*s1*s + g*s1*e2 + s1*s*e2 + s1*e2^2))/((e1 + g)*(e1 + s)*(g + e2)*(s + e2)*(n2 + n1))
  return(out)
}
betas = exp(fit.weather$par["c.0"]
            +fit.weather$par["c.p"]*weather.school$precipitation
            +fit.weather$par["c.t"]*weather.school$temperature
            +fit.weather$par["c.s"]*weather.school$school)
rts.base <- data.frame(week = weather.school$week,
                       x100 = rt(x100.base, fit.base$par[1]*(1+fit.base$par[2]*cos(2*pi*weather.school$week/52.14286+fit.base$par[3])),known)*c(rep(1,1055),rep(1-100/100,11),rep(1,58)),
                       x75 = rt(x75.base, fit.base$par[1]*(1+fit.base$par[2]*cos(2*pi*weather.school$week/52.14286+fit.base$par[3])),known)*c(rep(1,1055),rep(1-75/100,11),rep(1,58)),
                       x50 = rt(x50.base, fit.base$par[1]*(1+fit.base$par[2]*cos(2*pi*weather.school$week/52.14286+fit.base$par[3])),known)*c(rep(1,1055),rep(1-50/100,11),rep(1,58)))%>%
  pivot_longer(cols=!week,names_to="Efficacy",values_to="rt")%>%
  mutate(Efficacy=factor(Efficacy, levels=c("x50","x75","x100"), labels=c("50%","75%","100%")))
rts.weather <- data.frame(week = weather.school$week,
                          x100 = rt(x100.weather, betas*c(rep(1,1055),rep(1-100/100,11),rep(1,58)),known),
                          x75 = rt(x75.weather, betas*c(rep(1,1055),rep(1-75/100,11),rep(1,58)),known),
                          x50 = rt(x50.weather, betas*c(rep(1,1055),rep(1-50/100,11),rep(1,58)),known))%>%
  pivot_longer(cols=!week,names_to="Efficacy",values_to="rt")%>%
  mutate(Efficacy=factor(Efficacy, levels=c("x50","x75","x100"), labels=c("50%","75%","100%")))

weather.plot3 <- rts.weather %>% 
  mutate(year=as.Date(date_decimal(week/(365.25/7)+2000))) %>%
  filter(year>=dmy("1/1/2020"))%>%
  ggplot(aes(x=year, y=rt, color=Efficacy))+geom_line()+
  scale_x_date(date_breaks="3 months", date_labels="%b %Y", minor_breaks="1 month")+
  geom_vline(xintercept=dmy("23/3/2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("27/4/2020"), y=2.45, label="Lockdown", size=3)+
  geom_vline(xintercept=dmy("1.6.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.7.2020"), y=2.45, label="Winter", size=3)+
  geom_vline(xintercept=dmy("1.9.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.10.2020"), y=2.45, label="Spring", size=3)+
  geom_vline(xintercept=dmy("1.12.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.1.2021"), y=2.45, label="Summer", size=3)+
  geom_vline(xintercept=dmy("1.3.2021"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.4.2021"), y=2.45, label="Autumn", size=3)+
  geom_vline(xintercept=dmy("1.6.2021"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.7.2021"), y=2.45, label="Winter", size=3)+
  labs(x="", y=expression(R[eff]))+
  ggtitle("Weather model")+
  scale_color_brewer(palette="Dark2", name="Lockdown efficacy")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

base.plot3 <- rts.base %>% 
  mutate(year=as.Date(date_decimal(week/(365.25/7)+2000))) %>%
  filter(year>=dmy("1/1/2020"))%>%
  ggplot(aes(x=year, y=rt, color=Efficacy))+geom_line()+
  scale_x_date(date_breaks="3 months", date_labels="%b %Y", minor_breaks="1 month")+
  geom_vline(xintercept=dmy("23/3/2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("27/4/2020"), y=2, label="Lockdown", size=3)+
  geom_vline(xintercept=dmy("1.6.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.7.2020"), y=2, label="Winter", size=3)+
  geom_vline(xintercept=dmy("1.9.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.10.2020"), y=2, label="Spring", size=3)+
  geom_vline(xintercept=dmy("1.12.2020"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.1.2021"), y=2, label="Summer", size=3)+
  geom_vline(xintercept=dmy("1.3.2021"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.4.2021"), y=2, label="Autumn", size=3)+
  geom_vline(xintercept=dmy("1.6.2021"), linetype="dashed")+
  annotate(geom="text", x=dmy("15.7.2021"), y=2, label="Winter", size=3)+
  labs(x="", y=expression(R[eff]))+
  ggtitle("Sinusoidal model")+
  scale_color_brewer(palette="Dark2", name="Lockdown efficacy")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

arranged.plot3 <- ggarrange(base.plot3, weather.plot3, nrow=2, ncol=1, common.legend=TRUE, legend="bottom")
arranged.plot3
ggsave(filename="Figure6.tif",
       device="tiff",
       plot=arranged.plot3,
       path="plots",
       width=15,
       height=15,
       units="cm",
       dpi=300)
