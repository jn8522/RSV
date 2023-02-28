source("src/fit.R")
source("src/ODE.R")
source("src/plot.R")
library(ggpubr)
theme_set(theme_bw())
#reading in the model fits
fits.base <- readRDS("models/fits.base.CV.Rdata")
fits.weather <- readRDS("models/fits.weather.CV.Rdata")

#reading in some other relevant data
RSV <- read.csv("data/processed/RSV.csv")
known <- c(read.csv("data/raw/known.values.csv")$value, RSV$cases[1])
weather.school <- read.csv("data/processed/weather.school.csv")%>%filter(week<=nrow(RSV))
include=c(TRUE, #c.p
          TRUE, #c.t
          FALSE, #c.v
          TRUE) #c.s

#calculating the MSEs
MSE <- data.frame(i=1:5,
                  base=c(mse.parameters.base(fits.base[1,1]$par, data=RSV, known=known, i=1, k=5),
                         mse.parameters.base(fits.base[1,2]$par, data=RSV, known=known, i=2, k=5),
                         mse.parameters.base(fits.base[1,3]$par, data=RSV, known=known, i=3, k=5),
                         mse.parameters.base(fits.base[1,4]$par, data=RSV, known=known, i=4, k=5),
                         mse.parameters.base(fits.base[1,5]$par, data=RSV, known=known, i=5, k=5)),
                  weather=c(mse.parameters.weather(fits.weather[1,1]$par, data=RSV, known=known, weather=weather.school, include=include, i=1, k=5),
                            mse.parameters.weather(fits.weather[1,2]$par, data=RSV, known=known, weather=weather.school, include=include, i=2, k=5),
                            mse.parameters.weather(fits.weather[1,3]$par, data=RSV, known=known, weather=weather.school, include=include, i=3, k=5),
                            mse.parameters.weather(fits.weather[1,4]$par, data=RSV, known=known, weather=weather.school, include=include, i=4, k=5),
                            mse.parameters.weather(fits.weather[1,5]$par, data=RSV, known=known, weather=weather.school, include=include, i=5, k=5)))
avg.MSE <- c(base=mean(MSE$base), weather=mean(MSE$weather))

#generating figure 2 -- goodness of fit for each model
plots.base <- list()
legend.pos=c(0.9,0.9)
legend.margin=margin(t=1,r=3,b=1,l=1)
plot.margin=margin(t=0,r=0,b=0,l=5)
for (i in 1:5) {
  test.startyear <- (round((i-1)*1/5*nrow(RSV))+1)/52.14286+2000
  test.endyear <-(round(i*1/5*nrow(RSV)))/52.14286+2000
  plots.base[[i]] <- plot.fit.base(observed.df=RSV, variable=fits.base[1,i]$par, known=known)+
    geom_vline(xintercept=test.startyear,linetype="dashed")+
    geom_vline(xintercept=test.endyear,linetype="dashed")+
    theme(plot.margin=plot.margin)
  if(i==1){
    plots.base[[i]] <- plots.base[[i]]+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases), label="Test", size=3)+
      annotate("text", x=test.startyear+1.4, y=max(RSV$cases)-7,
               label=paste0("MSE=", round(MSE$base[i],1)), size=2)+
      annotate("text", x=test.endyear+0.75, y=max(RSV$cases), label="Train", size=3)+
      labs(title="Base Model",
           y="")+
      theme(plot.title = element_text(hjust = 0.5,
                                      margin=margin(t=5,b=5)),
            legend.position=legend.pos,
            legend.text=element_text(size=5),
            legend.margin=legend.margin)
  }
  else if(i==5){
    plots.base[[i]] <- plots.base[[i]]+
      annotate("text", x=test.startyear-0.75, y=max(RSV$cases), label="Train", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases), label="Test", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases)-7,
               label=paste0("MSE=", round(MSE$base[i],1)), size=2)+
      labs(caption=paste("Average MSE =",round(avg.MSE["base"],2)),
           y="")+
      theme(plot.caption = element_text(hjust = 0.5,
                                        margin=margin(t=-5,b=5)),
            legend.position="none")
  }else if(i==3){
    plots.base[[i]] <- plots.base[[i]]+
      annotate("text", x=test.startyear-0.75, y=max(RSV$cases), label="Train", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases), label="Test", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases)-7,
               label=paste0("MSE=", round(MSE$base[i],1)), size=2)+
      annotate("text", x=test.endyear+0.75, y=max(RSV$cases), label="Train", size=3)+
      theme(legend.position="none")
  }
  else{
    plots.base[[i]] <- plots.base[[i]]+
      annotate("text", x=test.startyear-0.75, y=max(RSV$cases), label="Train", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases), label="Test", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases)-7,
               label=paste0("MSE=", round(MSE$base[i],1)), size=2)+
      annotate("text", x=test.endyear+0.75, y=max(RSV$cases), label="Train", size=3)+
      labs(y="")+
      theme(legend.position="none")
  }
}
plots.weather <- list()
for (i in 1:5) {
  test.startyear <- (round((i-1)*1/5*nrow(RSV))+1)/52.14286+2000
  test.endyear <-(round(i*1/5*nrow(RSV)))/52.14286+2000
  plots.weather[[i]] <- plot.fit.weather(observed.df=RSV,
                                         variable=fits.weather[1,i]$par,
                                         known=known,
                                         weather=weather.school,
                                         include=include)+
    geom_vline(xintercept=test.startyear,linetype="dashed")+
    geom_vline(xintercept=test.endyear,linetype="dashed")+
    rremove("ylab")+
    theme(plot.margin=plot.margin,
          legend.position="none")
  if(i==1){
    plots.weather[[i]] <- plots.weather[[i]]+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases), label="Test", size=3)+
      annotate("text", x=test.startyear+1.45, y=max(RSV$cases)-7,
               label=paste0("MSE=", round(MSE$weather[i],1)), size=2)+
      annotate("text", x=test.endyear+0.75, y=max(RSV$cases), label="Train", size=3)+
      labs(title="Weather Model")+
      theme(plot.title = element_text(hjust = 0.5,
                                      margin=margin(t=5,b=5)),
            legend.position=legend.pos,
            legend.text=element_text(size=5),
            legend.margin=legend.margin)
  }
  else if(i==5){
    plots.weather[[i]] <- plots.weather[[i]]+
      annotate("text", x=test.startyear-0.75, y=max(RSV$cases), label="Train", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases), label="Test", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases)-7,
               label=paste0("MSE=", round(MSE$weather[i],1)), size=2)+
      labs(caption=paste("Average MSE =",round(avg.MSE["weather"],2)))+
      theme(plot.caption = element_text(hjust = 0.5,
                                        margin=margin(t=-5,b=5)),
            legend.position="none")
  }else{
    plots.weather[[i]] <- plots.weather[[i]]+
      annotate("text", x=test.startyear-0.75, y=max(RSV$cases), label="Train", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases), label="Test", size=3)+
      annotate("text", x=test.startyear+1.29452, y=max(RSV$cases)-7,
               label=paste0("MSE=", round(MSE$weather[i],1)), size=2)+
      annotate("text", x=test.endyear+0.75, y=max(RSV$cases), label="Train", size=3)+
      theme(legend.position="none")
  }
}

#putting them together
arranged <- ggarrange(plots.base[[1]], plots.weather[[1]],
                      plots.base[[2]], plots.weather[[2]],
                      plots.base[[3]], plots.weather[[3]],
                      plots.base[[4]], plots.weather[[4]],
                      plots.base[[5]], plots.weather[[5]],
                      heights=c(1,0.9,0.9,0.9,1),
                      nrow=5, ncol=2, 
                      common.legend=FALSE)
ggsave(filename="Figure2.tif",
       device="tiff",
       plot=arranged,
       path="plots",
       width=2250,
       height=2625,
       units="px",
       dpi=300)
