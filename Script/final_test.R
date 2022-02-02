#calculating the MSE of the base model on the test data
source("~/rsv-modelling/Script/Functions/ODE.R")
source("~/rsv-modelling/Script/Functions/training_preamble.R")

#usual preamble: training data, weather data, etc.
training_data()

#reading in test data
test.df <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/test.csv")
weather.train <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.train.csv")
weather.test <- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.test.csv")
weather <- union(weather.train, weather.test) %>%
  mutate(precipitation = precipitation/mean(weather.train$precipitation),
         evaporation = evaporation/mean(weather.train$evaporation),
         temperature = temperature/mean(weather.train$temperature),
         dewpoint = dewpoint/mean(weather.train$dewpoint),
         humidity = humidity/mean(weather.train$humidity))

#reading in the fitted parameters
par.base = read.csv("~/rsv-modelling/Output/Data/base.par.csv", header=FALSE)$V1
par.loglinear <- read.csv("~/rsv-modelling/Output/Data/Fits/Backwards selection/Loglinear round 1/Fit_loglinear_v_par.csv")$x
include=c(TRUE,TRUE,FALSE,TRUE)

#solving the ODE model
output.base <- ode.solve(variable=par.base,
                    known=known,
                    weeks = max(test.df$week),
                    noise=FALSE)
output.loglinear <- ode.solve.loglinear(variable=par.loglinear,
                              known=known,
                              noise=FALSE,
                              weather=weather,
                              include=include)

#calculating mse on test data
mse.base.df <- inner_join(test.df, output.base, by="week")
mse.base <- mean((mse.base.df$cases.x-mse.base.df$cases.y)^2)
write.csv(mse.base, file="~/rsv-modelling/Output/Data/test_mse_base.csv")

mse.loglinear.df <- inner_join(test.df, output.loglinear, by="week")
mse.loglinear <- mean((mse.loglinear.df$cases.x-mse.loglinear.df$cases.y)^2)
write.csv(mse.loglinear, file="~/rsv-modelling/Output/Data/test_mse_extended.csv")

#plotting fit over all the data
plot.base.df <- union(train.df, test.df, by="week")%>%
  left_join(output.base, by="week")%>%
  rename(Observed=cases.x, Modelled=cases.y)%>%
  pivot_longer(cols=!week, names_to="source", values_to="cases")%>%
  mutate(year=week/(365.25/7)+2000)

grey="grey30"
base.plot <- plot.base.df %>% 
  ggplot(aes(x=year, y=cases, color=source))+
  geom_line()+
  scale_x_continuous(breaks=seq(0,ceiling(max(plot.base.df$year)),by=1))+
  scale_color_manual(values=c("red", grey))+
  theme(legend.title=element_blank(), legend.position="none",
        plot.title=element_text(size=10), axis.title.y = element_text(size=9))+
  ggtitle("Sinusoidal model")+
  labs(x="", y="Cases")+
  geom_vline(xintercept=test.df$week[1]/(365.25/7)+2000, linetype="dashed")+
  annotate("text", x=(test.df$week[1]-30)/(365.25/7)+2000, y=max(plot.base.df$cases), label="Train", size=3)+
  annotate("text", x=(test.df$week[1]+30)/(365.25/7)+2000, y=max(plot.base.df$cases), label="Test", size=3)+
  annotate("text", x=(test.df$week[1]+63)/(365.25/7)+2000, y=max(plot.base.df$cases)-5, label=paste("MSE = ", round(mse.base,2)), size=3)

plot.loglinear.df <- union(train.df, test.df, by="week")%>%
  left_join(output.loglinear, by="week")%>%
  rename(Observed=cases.x, Modelled=cases.y)%>%
  pivot_longer(cols=!week, names_to="source", values_to="cases")%>%
  mutate(year=week/(365.25/7)+2000)

loglinear.plot <- plot.loglinear.df %>% 
  ggplot(aes(x=year, y=cases, color=source))+
  geom_line()+
  scale_x_continuous(breaks=seq(0,ceiling(max(plot.loglinear.df$year)),by=1))+
  scale_color_manual(values=c("red", grey))+
  theme(legend.title=element_blank(), legend.position="none",
        plot.title=element_text(size=10), axis.title.y = element_text(size=9))+
  ggtitle("Extended model")+
  labs(x="", y="Cases")+
  geom_vline(xintercept=test.df$week[1]/(365.25/7)+2000, linetype="dashed")+
  annotate("text", x=(test.df$week[1]-30)/(365.25/7)+2000, y=max(plot.loglinear.df$cases), label="Train", size=3)+
  annotate("text", x=(test.df$week[1]+30)/(365.25/7)+2000, y=max(plot.loglinear.df$cases), label="Test", size=3)+
  annotate("text", x=(test.df$week[1]+63)/(365.25/7)+2000, y=max(plot.loglinear.df$cases)-5, label=paste("MSE = ", round(mse.loglinear,2)), size=3)

ggsave(filename="base_evaluate.png",
       plot=base.plot,
       path="~/rsv-modelling/Output/Plots",
       width=19,
       height=10.7,
       units="cm",
       dpi=300)

ggsave(filename="loglinear_evaluate.png",
       plot=loglinear.plot,
       path="~/rsv-modelling/Output/Plots",
       width=19,
       height=10.7,
       units="cm",
       dpi=300)

library(ggpubr)
arranged <- ggarrange(base.plot, loglinear.plot, nrow=2, ncol=1, common.legend = TRUE, legend="bottom")
ggsave(filename="evaluate_both.tif",
       device="tiff",
       plot=arranged,
       path="~/rsv-modelling/Output/Plots",
       width=15,
       height=15,
       units="cm",
       dpi=300)

