source("src/fit.R")
source("src/ODE.R")
source("src/plot.R")
library(ggpubr)
theme_set(theme_bw())
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

#functions that calculate the reproduction numbers
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

#solving base model ODE to get numbers of susceptibles
#determining initial conditions based on provided ratios
n1 = known[8]
n2 = known[9]
ra1 = fit.base$par[6]
ra2 = fit.base$par[7]
ra3 = fit.base$par[8]
ra4 = fit.base$par[9]
ra5 = fit.base$par[10]
ra6 = fit.base$par[11]
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
p=fit.base$par[5]
y0 = c(s1,e1,i1,r1,s2,e2,i2,r2,
       known[10]/p)

ode.output.base <- ode(func=seir.model.base, 
                       y=y0, 
                       times=1:nrow(RSV),
                       parms=fit.base$par,
                       known=known) %>%
  as.data.frame()
names(ode.output.base) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new")

#and for weather model
n1 = known[8]
n2 = known[9]
ra1 = fit.weather$par[length(fit.weather$par)-5]
ra2 = fit.weather$par[length(fit.weather$par)-4]
ra3 = fit.weather$par[length(fit.weather$par)-3]
ra4 = fit.weather$par[length(fit.weather$par)-2]
ra5 = fit.weather$par[length(fit.weather$par)-1]
ra6 = fit.weather$par[length(fit.weather$par)]
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
p = fit.base$par[2]
y0 = c(s1,e1,i1,r1,s2,e2,i2,r2,
       known[10]/p)
cs <- rep(0, 4)
c.0 <- fit.weather$par[3]
index=0
for(i in 1:4){
  if(include[i]==TRUE){
    cs[i] = fit.weather$par[i+3-index]
  }
  else{
    index = index+1
  }
}
betas = exp(c.0
            +cs[1]*weather.school$precipitation
            +cs[2]*weather.school$temperature
            +cs[3]*weather.school$vapour
            +cs[4]*weather.school$school)

ode.output.ext <- ode(func=seir.model.weather, 
                      y=y0, 
                      times=1:nrow(RSV),
                      parms=fit.weather$par,
                      known=known,
                      betas=betas) %>%
  as.data.frame()

names(ode.output.ext) <- c("week", "S1", "E1", "I1", "R1", "S2", "E2", "I2", "R2", "new")

#plotting R0s
r0s <- data.frame(year = RSV$week/(365.25/7)+2000,
                  Sinusoidal = r0(ode.output.base, fit.base$par[1]*(1+fit.base$par[2]*cos(2*pi*RSV$week/52.14286+fit.base$par[3])),known),
                  Extended = r0(ode.output.ext, betas,known))

r0.plot <- r0s %>%
  pivot_longer(!year, names_to="Model", values_to = "R0")%>%
  ggplot(aes(x=year, y=R0, color=Model))+
  geom_line()+
  scale_x_continuous(breaks=seq(2000,2013,by=1))+
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0))+
  scale_color_brewer(palette="Dark2",
                     labels=c("Weather", "Sinusoidal"),
                     name="Model")+
  geom_hline(yintercept=1, linetype="dashed")+
  ggtitle("Basic reproduction number")+
  labs(x="",
       y=expression(R[0]))+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

#plotting effective reproduction number rt
rts <- data.frame(year = RSV$week/(365.25/7)+2000,
                  sinusoidal = rt(ode.output.base, fit.base$par[1]*(1+fit.base$par[2]*cos(2*pi*RSV$week/52.14286+fit.base$par[3])),known),
                  ext = rt(ode.output.ext, betas,known))

rt.plot <- rts %>%
  pivot_longer(cols=!year, names_to="variable", values_to="value")%>%
  ggplot(aes(x=year, y=value, color=variable))+
  geom_line()+
  scale_x_continuous(breaks=seq(2000,2013,by=1))+
  theme(legend.position = "none")+
  scale_color_brewer(palette="Dark2",
                     labels=c("Weather", "Sinusoidal"),
                     name="Model")+
  geom_hline(yintercept=1, linetype="dashed")+
  labs(x="",y=expression(R[eff]))+
  ggtitle("Effective reproduction number")+
  theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
        legend.title=element_text(size=8))

library(ggpubr)
arranged <- ggarrange(rt.plot, r0.plot, nrow=2, ncol=1, common.legend=TRUE, legend="bottom")
ggsave(filename="Figure3.tif",
       device="tiff",
       plot=arranged,
       path="plots",
       width=15,
       height=15,
       units="cm",
       dpi=300)
