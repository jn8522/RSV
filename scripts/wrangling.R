library(tidyverse)
library(lubridate)
library(purrr)
library(zoo)
#RSV DATA
#importing the raw data
flu_rsv <- read.csv("data/raw/Flu_RSV.csv")
RSV <- filter(flu_rsv, (rsv=="Yes"|rsvflu==1) & metrure=="Metro")%>%
  mutate(weeknum = interval(dmy("1.1.2000"),dmy(date)) %/% weeks(1)+1)%>%
  count(weeknum)

RSV.df <- data.frame(weeknum = 1:max(RSV$weeknum),
                       n = 0) %>%
  bind_rows(RSV) %>%
  group_by(weeknum) %>%
  summarise(cases = sum(n)) %>%
  rename(week=weeknum)

write_csv(RSV.df, "data/processed/RSV.csv")

#SCHOOL AND WEATHER DATA
#reading in the school data
school.data <- read.csv("data/raw/school_terms.csv") %>% 
  data.table::transpose(keep.names = "Year",
                        make.names = "Year") %>%
  mutate(Year = substring(Year, 2))

#wrangling them into a date format
school <- map(school.data[,2:15], ~dmy(paste0(.x,school.data[,1]))) %>%
  as.data.frame() %>%
  mutate(Year = dmy(paste0("1/1/",school.data$Year)))

#working out what week of the year they tend to fall on
interval(school$Year, school$start.term2) %/% weeks(1)+1
school.weeks <- map(school[,1:14], ~interval(school[,15],.x)%/%weeks(1))%>%
  as.data.frame()
dates <- school.weeks %>% summarise(across(everything(), list(min)))

#writing a function that takes a vector of week numbers and returns a vector of whether school is on or off
term <- function(weeks, dates) {
  out <- rep(0,length(weeks))
  for(i in 1:length(weeks)){
    if(dates$start.term1_1<=weeks[i]&
       weeks[i]<=dates$end.term1_1){
      out[i] <- 1
    }
    if(dates$start.term2_1<=weeks[i]&
       weeks[i]<=dates$end.term2_1){
      out[i] <- 1
    }
    if(dates$start.term3_1<=weeks[i]&
       weeks[i]<=dates$end.term3_1){
      out[i] <- 1
    }
    if(dates$start.term4_1<=weeks[i]&
       weeks[i]<=dates$end.term4_1){
      out[i] <- 1
    }
  }
  return(out)
}

#reading in the weather data
weather.data <- read.csv("data/raw/bom_data.txt")

#applying the school terms function to the study period
school.terms <- data.frame(date = dmy(weather.data$Day.Month.Year.in.DD.MM.YYYY.format)) %>%
  mutate(year.start = dmy(paste0("1/1/",year(date))),
         week.overall = interval(dmy("1.1.2000"),date) %/% weeks(1)+1) %>%
  mutate(week.year = interval(year.start,date) %/% weeks(1)+1) %>%
  mutate(school = term(week.year, dates))%>%
  group_by(week.overall)%>%
  summarise(school=round(mean(school))) %>%
  rename(week=week.overall)

#tidying the weather data
weather.data <- read.csv("data/raw/bom_data.txt") %>%
  select(date = Day.Month.Year.in.DD.MM.YYYY.format,
         precipitation = Precipitation.in.the.24.hours.before.9am..local.time..in.mm,
         precipitation.quality = Quality.of.precipitation.value,
         precipitation.acc = Accumulated.number.of.days.over.which.the.precipitation.was.measured,
         evaporation = Evaporation.in.24.hours.before.9am..local.time..in.mm,
         evaporation.quality = Quality.of.evaporation.in.24.hours.before.9am..local.time.,
         evaporation.acc = Days.of.accumulation.for.evaporation,
         temperature = Air.temperature.observation.at.09.hours.Local.Time.in.Degrees.C,
         temperature.quality = Quality.of.air.temperature.observation.at.09.hours.Local.Time,
         dewpoint = Dew.point.temperature.observation.at.09.hours.Local.Time.in.Degrees.C,
         dewpoint.quality = Quality.of.dew.point.temperature.observation.at.09.hours.Local.Time,
         humidity = Relative.humidity.for.observation.at.09.hours.Local.Time.in.percentage..,
         humidity.quality = Quality.of.relative.humidity.for.observation.at.09.hours.Local.Time,
         date = Day.Month.Year.in.DD.MM.YYYY.format)

#Will use last observation carried forward for readings that are 'suspect' or 'wrong'
for (i in 1:nrow(weather.data)) {
  if(weather.data$precipitation.quality[i] %in% c("S", "W")){
    weather.data$precipitation[i] <- NA
  }
  if(weather.data$evaporation.quality[i] %in% c("S", "W")){
    weather.data$evaporation[i] <- NA
  }
  if(weather.data$temperature.quality[i] %in% c("S", "W")){
    weather.data$temperature[i] <- NA
  }
  if(weather.data$dewpoint.quality[i] %in% c("S", "W")){
    weather.data$dewpoint[i] <- NA
  }
  if(weather.data$humidity.quality[i] %in% c("S", "W")){
    weather.data$humidity[i] <- NA
  }
  #de-accumulating evaporation readings that went over 2 days
  if(!is.na(weather.data$evaporation.acc[i])){
    if(weather.data$evaporation.acc[i]==2){
      weather.data$evaporation[i] <- weather.data$evaporation[i]/2
      weather.data$evaporation[i-1] <- weather.data$evaporation[i]
    }
  }
}

#simplifying output
weather.data <- weather.data %>%
  select(-precipitation.acc, -evaporation.acc) %>%
  #last observation carried forward
  zoo::na.locf() %>%
  #grouping by week
  mutate(week = interval(dmy("1.1.2000"),dmy(date)) %/% weeks(1)+1,
         vapour = exp(1.8096+(17.269425*dewpoint)/(237.3+dewpoint))) %>%
  group_by(week) %>%
  #converting to weekly totals
  summarise(precipitation=sum(precipitation),
            evaporation=sum(evaporation),
            temperature=sum(temperature),
            dewpoint=sum(dewpoint),
            humidity=sum(humidity),
            vapour=sum(vapour))%>%
  #and normalising by the means
  mutate(precipitation = (precipitation-mean(precipitation))/sd(precipitation),
         evaporation = (evaporation-mean(evaporation)/sd(evaporation)),
         temperature = (temperature-mean(temperature))/sd(temperature),
         dewpoint = (dewpoint-mean(dewpoint))/sd(dewpoint),
         humidity = (humidity-mean(humidity))/sd(humidity),
         vapour=(vapour-mean(vapour))/sd(vapour))
#joining this with the school data
weather.school <- left_join(weather.data, school.terms, by="week")

#outputting the results
write_csv(weather.school, "data/processed/weather.school.csv")