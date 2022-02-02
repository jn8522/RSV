#usual preamble when fitting models to training data
training_data <- function() {
  #reading in training data
  train.df <<- read.csv("~/rsv-modelling/Output/Data/Wrangled data/train.csv")
  
  #setting these as known values for the SEIR model
  known.values = read.csv("~/rsv-modelling/Data/known.values.csv")
  known <<- c(known.values$value, train.df$cases[1])
  
  #reading in the weather data
  weather <<- read.csv("~/rsv-modelling/Output/Data/Wrangled data/weather.train.csv") %>%
    mutate(precipitation = precipitation/mean(precipitation),
           evaporation = evaporation/mean(evaporation),
           temperature = temperature/mean(temperature),
           dewpoint = dewpoint/mean(dewpoint),
           humidity = humidity/mean(humidity))
  
}