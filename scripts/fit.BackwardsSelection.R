#Fitting the weather model to the training data
library(dfoptim)
library(parallel)

source("src/fit.R")
source("src/ODE.R")
source("src/plot.R")

#reading in data
RSV <- read.csv("data/processed/RSV.csv")
known <- c(read.csv("data/raw/known.values.csv")$value, RSV$cases[1])
weather.school <- read.csv("data/processed/weather.school.csv")%>%filter(week<=nrow(RSV))

#choosing initial values
initial = c(nu=0.03,
            p=0.01,
            c.0=1.4,
            c.p=0.08,
            c.t=-0.3,
            c.v=0.02,
            c.s=0.1,
            ra1=0.01,
            ra2=0,
            ra3=0.3,
            ra4=0,
            ra5=0.3,
            ra6=0.05)

#setting the lower and upper bounds for coefficients c.p, c.t, c.v and c.s
lower.bounds=c(nu=0, 
        p=0,
        c.0=-Inf,
        c.p=-Inf,
        c.t=-Inf,
        c.v=-Inf,
        c.s=-Inf,
        ra1=0,ra2=0,ra3=0,ra4=0,ra5=0,ra6=0)
upper.bounds=c(nu=Inf,
        p=1,
        c.0=Inf,
        c.p=Inf,
        c.t=Inf,
        c.v=Inf,
        c.s=Inf,
        ra1=1,ra2=1,ra3=1,ra4=1,ra5=1,ra6=1)
include=list(c(TRUE,TRUE,TRUE,TRUE), #all variables
             c(FALSE,TRUE,TRUE,TRUE), #dropping precipitation
             c(TRUE,FALSE,TRUE,TRUE), #dropping temperature
             c(TRUE,TRUE,FALSE,TRUE), #dropping vapour pressure
             c(TRUE,TRUE,TRUE,FALSE)) #dropping school terms
fits.backwards.selection <- mcmapply(hjkb,
                                     mc.cores = 5,
                                     include=include,
                                     par=list(initial,
                                              initial[-4],
                                              initial[-5],
                                              initial[-6],
                                              initial[-7]),
                                     lower=list(lower.bounds,
                                                lower.bounds[-4],
                                                lower.bounds[-5],
                                                lower.bounds[-6],
                                                lower.bounds[-7]),
                                     upper=list(upper.bounds,
                                                upper.bounds[-4],
                                                upper.bounds[-5],
                                                upper.bounds[-6],
                                                upper.bounds[-7]),
                                     MoreArgs=list(fn=neglk.parameters.weather,
                                                   control=list(info=TRUE),
                                                   data=RSV,
                                                   known=known,
                                                   weather=weather.school,
                                                   i=0,
                                                   k=0)) #setting k=0 so we're fitting to the entire time series (not using CV)

#fitting one more time to check they all *really* converged
fits.backwards.selection2 <- mcmapply(hjkb,
                                      mc.cores = 5,
                                      include=include,
                                      par=list(fits.backwards.selection[1,1]$par,
                                               fits.backwards.selection[1,2]$par,
                                               fits.backwards.selection[1,3]$par,
                                               fits.backwards.selection[1,4]$par,
                                               fits.backwards.selection[1,5]$par),
                                      lower=list(lower.bounds,
                                                 lower.bounds[-4],
                                                 lower.bounds[-5],
                                                 lower.bounds[-6],
                                                 lower.bounds[-7]),
                                      upper=list(upper.bounds,
                                                 upper.bounds[-4],
                                                 upper.bounds[-5],
                                                 upper.bounds[-6],
                                                 upper.bounds[-7]),
                                      MoreArgs=list(fn=neglk.parameters.weather,
                                                    data=RSV,
                                                    known=known,
                                                    weather=weather.school,
                                                    i=0, 
                                                    k=0))%>% #setting k=0 so we're fitting to the entire time series (not using CV)
  as.data.frame()
colnames(fits.backwards.selection2) <- c("All", "-P", "-T", "-V", "-S")

#calculating AICs
L <- c(All=fits.backwards.selection2["value",1]$value,
       `-P`=fits.backwards.selection2["value",2]$value,
       `-T`=fits.backwards.selection2["value",3]$value,
       `-V`=fits.backwards.selection2["value",4]$value,
       `-S`=fits.backwards.selection2["value",5]$value)
k <- c(fits.backwards.selection2["par",1]$par %>% length(),
       fits.backwards.selection2["par",2]$par %>% length(),
       fits.backwards.selection2["par",3]$par %>% length(),
       fits.backwards.selection2["par",4]$par %>% length(),
       fits.backwards.selection2["par",5]$par %>% length())
AIC <- 2*k+2*t(L)
fits.backwards.selection2["AIC",] <- AIC

#doing a second round dropping vapour pressure and one other variable
fits.backwards.selection3 <- mcmapply(hjkb,
                                      mc.cores = 3,
                                      include=list(c(FALSE,TRUE,FALSE,TRUE), #dropping precipitation
                                                   c(TRUE,FALSE,FALSE,TRUE), #dropping temperature
                                                   c(TRUE,TRUE,FALSE,FALSE)), #dropping school terms
                                      par=list(`-V,-P`=fits.backwards.selection2$`-V`$par[-4],
                                               `-V,-T`=fits.backwards.selection2$`-V`$par[-5],
                                               `-V,-S`=fits.backwards.selection2$`-V`$par[-6]),
                                      lower=list(lower.bounds[-6][-4],
                                                 lower.bounds[-6][-5],
                                                 lower.bounds[-6][-6]),
                                      upper=list(upper.bounds[-6][-4],
                                                 upper.bounds[-6][-5],
                                                 upper.bounds[-6][-6]),
                                      MoreArgs=list(fn=neglk.parameters.weather,
                                                    data=RSV,
                                                    known=known,
                                                    weather=weather.school,
                                                    i=0, 
                                                    k=0))%>% #setting k=0 so we're fitting to the entire time series (not using CV)
  as.data.frame()

#checking that they definitely converged
fits.backwards.selection4 <- mcmapply(hjkb,
                                      mc.cores = 3,
                                      include=list(c(FALSE,TRUE,FALSE,TRUE), #dropping precipitation
                                                   c(TRUE,FALSE,FALSE,TRUE), #dropping temperature
                                                   c(TRUE,TRUE,FALSE,FALSE)), #dropping school terms
                                      par=list(fits.backwards.selection3$V1$par,
                                               fits.backwards.selection3$V2$par,
                                               fits.backwards.selection3$V3$par),
                                      lower=list(lower.bounds[-6][-4],
                                                 lower.bounds[-6][-5],
                                                 lower.bounds[-6][-6]),
                                      upper=list(upper.bounds[-6][-4],
                                                 upper.bounds[-6][-5],
                                                 upper.bounds[-6][-6]),
                                      MoreArgs=list(fn=neglk.parameters.weather,
                                                    data=RSV,
                                                    known=known,
                                                    weather=weather.school,
                                                    i=0, 
                                                    k=0))%>% #setting k=0 so we're fitting to the entire time series (not using CV)
  as.data.frame()
colnames(fits.backwards.selection4) <- c("-V,P", "-V,T", "-V,S")

#calculating AICs
L <- c(`-V,P`=fits.backwards.selection4["value",1]$value,
       `-V,T`=fits.backwards.selection4["value",2]$value,
       `-V,S`=fits.backwards.selection4["value",3]$value)
k <- c(fits.backwards.selection4["par",1]$par %>% length(),
       fits.backwards.selection4["par",2]$par %>% length(),
       fits.backwards.selection4["par",3]$par %>% length())
AIC <- 2*k+2*t(L)
fits.backwards.selection4["AIC",] <- AIC
fits.backwards.selection.all <- bind_cols(fits.backwards.selection2,fits.backwards.selection4)
saveRDS(fits.backwards.selection.all, file="models/fits.backwards.selection.Rdata")
