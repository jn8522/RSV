#overlay observed data with model output for parameter values theta
#for the base model:
plot.fit.base <- function(observed.df, #dataframe of observed data
                     variable, #variable parameters
                     known #fixed parameters
) {
  fitted.df <- ode.solve.base(variable=variable,
                         known=known,
                         weeks=nrow(observed.df),
                         noise=FALSE)
  combined.df <- left_join(observed.df, fitted.df, 
                           by="week",
                           suffix=c(".observed", ".fitted"))%>%
    rename(Observed=cases.observed, Modelled=cases.fitted)%>%
    pivot_longer(!week, names_to = "source", values_to="value")%>%
    mutate(year = week/52.14286+2000)
  end_year = max(combined.df$year) %>% ceiling()
  plot <- combined.df %>%
    ggplot(aes(x=year, y=value, color=source))+geom_line()+
    scale_x_continuous(breaks=seq(0,end_year,by=2))+
    scale_color_manual(values=c("red", "#00000080"))+
    theme(legend.title=element_blank(), legend.position="bottom")+
    labs(x="", y="Weekly hospitalizations")+
    theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
          legend.title=element_blank())
  return(plot)
}

#for the extended model
plot.fit.weather <- function(observed.df, #dataframe of observed data
                          variable, #variable parameters
                          known, #fixed parameters
                          weather, #weather and school data
                          include #which variables to include
) {
  fitted.df <- ode.solve.weather(variable=variable,
                              known=known,
                              noise=FALSE,
                              weather=weather,
                              include=include)
  combined.df <- left_join(observed.df, fitted.df, 
                           by="week",
                           suffix=c(".observed", ".fitted"))%>%
    rename(Observed=cases.observed, Modelled=cases.fitted)%>%
    pivot_longer(!week, names_to = "source", values_to="value")%>%
    mutate(year = week/52.14286+2000)
  end_year = max(combined.df$year) %>% ceiling()
  plot <- combined.df %>%
    ggplot(aes(x=year, y=value, color=source))+geom_line()+
    scale_x_continuous(breaks=seq(0,end_year,by=2))+
    scale_color_manual(values=c("red", "#00000080"))+
    theme(legend.title=element_blank(), legend.position="bottom")+
    labs(x="", y="Weekly hospitalizations")+
    theme(plot.title=element_text(size=10), axis.title.y = element_text(size=9),
          legend.title=element_blank())
  return(plot)
}