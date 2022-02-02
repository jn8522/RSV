library(tidyverse)
library(lubridate)
library(purrr)

#reading in term dates
data <- read.csv("~/rsv-modelling/Data/school_terms.csv") %>% 
  data.table::transpose(keep.names = "Year",
                        make.names = "Year") %>%
  mutate(Year = substring(Year, 2))

#wrangling them into a date format
school <- map(data[,2:15], ~dmy(paste0(.x,data[,1]))) %>%
  as.data.frame() %>%
  mutate(Year = dmy(paste0("1/1/",data$Year)))

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

#applying this to the study period
#reading in the weather data
weather <- read.csv("~/rsv-modelling/Data/bom_data.txt")

school.terms <- data.frame(date = dmy(weather$Day.Month.Year.in.DD.MM.YYYY.format)) %>%
  mutate(year.start = dmy(paste0("1/1/",year(date))),
         week.overall = interval(dmy("1.1.2000"),date) %/% weeks(1)+1) %>%
  mutate(week.year = interval(year.start,date) %/% weeks(1)+1) %>%
  mutate(school = term(week.year, dates))%>%
  group_by(week.overall)%>%
  summarise(school=round(mean(school))) %>%
  rename(week=week.overall)

write_csv(school.terms, "~/rsv-modelling/Output/Data/Wrangled data/school.csv")