#plotting fits
#usual preamble
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")
training_data()

#base model
par.base <- read.csv("~/rsv-modelling/Output/Data/base.par.csv", header=FALSE)$V1

fit.plot.base <- plot.fit(train.df,
                         variable=par.base,
                         known=known)+ggtitle("Base model")

ggsave(filename="Fit_base_actual.png",
       plot=fit.plot.base,
       path="~/rsv-modelling/Output/Plots",
       width=24,
       height=13.5,
       units="cm",
       dpi=300)

mse.base <- ssq.parameters(par.base, train.df, known)
aic.base <- 2*length(par.base)+2*neglk.parameters(par.base, train.df,known)
  
#simplified loglinear model
par.loglinear <- read.csv("~/rsv-modelling/Output/Data/Fits/Backwards selection/Loglinear round 1/Fit_loglinear_v_par.csv")$x
include=c(TRUE,TRUE,FALSE,TRUE)

fit.plot.loglinear <- plot.fit.loglinear(train.df,
                                         variable=par.loglinear,
                                         known=known,
                                         weather=weather,
                                         include=include)+ggtitle("Extended model")
fit.plot.loglinear
ggsave(filename="Fit_loglinear_simplified.png",
       plot=fit.plot.loglinear,
       path="~/rsv-modelling/Output/Plots",
       width=24,
       height=13.5,
       units="cm",
       dpi=300)

#plotting transmissibility
library(zoo)
library(RColorBrewer)

weather.train <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.train.csv")
weather.test <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.test.csv")
weather <- union(weather.train, weather.test) %>%
  mutate(precipitation = precipitation/mean(weather.train$precipitation),
         evaporation = evaporation/mean(weather.train$evaporation),
         temperature = temperature/mean(weather.train$temperature),
         dewpoint = dewpoint/mean(weather.train$dewpoint),
         humidity = humidity/mean(weather.train$humidity))


k=10
transmissibility <- data.frame(year=weather$week/(365.25/7),
                               `Extended (10-day rolling average)` = rollmean(exp(par.loglinear[3]
                                               +par.loglinear[4]*weather$precipitation
                                               +par.loglinear[5]*weather$temperature
                                               +par.loglinear[6]*weather$school),k=k,fill=NA),
                              Sinusoidal=par.base[1]*(1+par.base[2]*cos(2*pi*weather$week/52.14286+par.base[3])))
beta.plot <- transmissibility %>% 
  pivot_longer(!year, names_to="Model", values_to = "Beta")%>%
  ggplot(aes(x=year, y=Beta, color=Model))+
  geom_line()+
  labs(x="Year", y="Beta")+
  scale_x_continuous(breaks=seq(0,13,by=1))+
  theme(legend.position = "bottom")+
  scale_color_brewer(palette="Dark2")
beta.plot

ggsave(filename="transmissibility.png",
       plot=beta.plot,
       path="~/rsv-modelling/Output/Plots",
       width=24,
       height=13.5,
       units="cm",
       dpi=300)

#school vs other things
school.plot <- data.frame(cases=train.df$cases,
           school=weather.train$school,
           year=train.df$week/(365.25/7),
           precipitation=weather.train$precipitation,
           temperature=weather.train$temperature) %>%
  filter(year<=3)%>%
  pivot_longer(!year, names_to="variable",values_to="value") %>%
  ggplot(aes(x=year,y=value))+
  geom_line()+
  facet_grid(rows=vars(variable),scales="free_y")
school.plot
ggsave(filename="school.png",
       plot=school.plot,
       path="~/rsv-modelling/Output/Plots",
       width=24,
       height=13.5,
       units="cm",
       dpi=300)

#looking at relative contribution of each variable
predictors <- data.frame(precipitation=weather.train$precipitation,
                   temperature=weather.train$temperature,
                   school=weather.train$school)
beta=exp(par.loglinear[3]
         +par.loglinear[4]*weather.train$precipitation
         +par.loglinear[5]*weather.train$temperature
         +par.loglinear[6]*weather.train$school)
library(hier.part)
hier <- hier.part(beta,predictors,family="poisson",gof="RMSPE",barplot=FALSE)

#reproduction numbers
r0 <- function(data, beta){
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

rt <- function(data, beta){
  s1 = data$S1
  s2 = data$S2
  n1=data$S1+data$E1+data$I1+data$R1
  n2=data$S2+data$E2+data$I2+data
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

#solving base model ODE to get numbers of susceptibles
#determining initial conditions based on provided ratios
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
       known[10]/p)

ode.output.base <- ode(func=seir.model, 
                  y=y0, 
                  times=1:nrow(weather),
                  parms=par.base,
                  known=known) %>%
  as.data.frame()
names(ode.output.base) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new")

#and for extended model
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
       known[10]/p)
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

names(ode.output.ext) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new")

#plotting R0s
r0s <- data.frame(year = weather$week/(365.25/7)+2000,
                  Sinusoidal = r0(ode.output.base, par.base[1]*(1+par.base[2]*cos(2*pi*weather$week/52.14286+par.base[3]))),
                  Extended = r0(ode.output.ext, betas))

r0.plot <- r0s %>%
  pivot_longer(!year, names_to="Model", values_to = "R0")%>%
  ggplot(aes(x=year, y=R0, color=Model))+
  geom_line()+
  scale_x_continuous(breaks=seq(2000,2013,by=1))+
  theme(legend.position = "none",
        plot.caption = element_text(hjust = 0),
        plot.title=element_text(size=10),
        axis.title.y = element_text(size=9))+
  scale_color_brewer(palette="Dark2")+
  geom_hline(yintercept=1, linetype="dashed")+
  ggtitle("Basic reproduction number")+
  labs(x="",
       y=expression(R[0]))
r0.plot

ggsave(filename="r0.png",
       plot=r0.plot,
       path="~/rsv-modelling/Output/Plots",
       width=24,
       height=13.5,
       units="cm",
       dpi=300)


#plotting effective reproduction number rt
rt <- function(data, beta){
  s1 = data$S1
  s2 = data$S2
  n1=data$S1+data$E1+data$I1+data$R1
  n2=data$S2+data$E2+data$I2+data
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

#solving base model ODE to get numbers of susceptibles
#determining initial conditions based on provided ratios
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
       known[10]/p)

#solving ODE and doing some data wrangling
ode.output.base <- ode(func=seir.model, 
                  y=y0, 
                  times=1:nrow(weather),
                  parms=par.base,
                  known=known) %>%
  as.data.frame()
names(ode.output.base) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new")

#and for extended model
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
       known[10]/p)
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

#solving ODE and doing some data wrangling
ode.output.ext <- ode(func=seir.model.ext, 
                  y=y0, 
                  times=1:nrow(weather),
                  parms=par.loglinear,
                  known=known,
                  betas=betas) %>%
  as.data.frame()

names(ode.output.ext) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new")

rts <- data.frame(year = weather$week/(365.25/7)+2000,
                  Sinusoidal = rt(ode.output.base, par.base[1]*(1+par.base[2]*cos(2*pi*weather$week/52.14286+par.base[3]))),
                  Extended = rt(ode.output.ext, betas))

rt.plot <- rts %>%
  pivot_longer(cols=!year, names_to="Model", values_to="value")%>%
  ggplot(aes(x=year, y=value, color=Model))+
  geom_line()+
  scale_x_continuous(breaks=seq(2000,2013,by=1))+
  theme(legend.position = "none",
        plot.title=element_text(size=10),
        axis.title.y = element_text(size=9))+
  scale_color_brewer(palette="Dark2")+
  geom_hline(yintercept=1, linetype="dashed")+
  labs(x="",y=expression(R[eff]))+
  ggtitle("Effective reproduction number")
rt.plot

ggsave(filename="rt.png",
       plot=rt.plot,
       path="~/rsv-modelling/Output/Plots",
       width=24,
       height=13.5,
       units="cm",
       dpi=300)

library(ggpubr)
arranged <- ggarrange(rt.plot, r0.plot, nrow=2, ncol=1, common.legend=TRUE, legend="bottom")
ggsave(filename="r0_rt.tif",
       device="tiff",
       plot=arranged,
       path="~/rsv-modelling/Output/Plots",
       width=15,
       height=15,
       units="cm",
       dpi=300)