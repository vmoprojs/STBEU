\name{OntarioTemp}
\alias{OntarioTemp}
\title{Ontario Hourly temperature}
\usage{data(OntarioSummer2018Hour)}
\description{
  Ontario Hourly temperature on 2018.
  
  
  
  
  # ***** (Not necesasary to run) Final Dataset used in STBEU was generated using the following:

# library(weathercan)
# library(dplyr)
# library(GeoModels)
# rm( list = ls())
# graphics.off()
# gc()
# cat("\014")
# 
# 
# mb <- filter(stations,
#              prov == "ON",
#              interval == "day",
#              end >= 2018) %>%
#   select(-prov)
# 
# 
# mb_weather_all <- weather_dl(station_ids = mb$station_id,
#                              start = "2018-03-20",
#                              interval = "hour", quiet = TRUE)
# 
# # Summer - June 21 to September 21 (YES)
# st <- "2018-06-21"
# ed <- "2018-09-21"
# 
# 
# fil <- mb_weather_all$date >= st & mb_weather_all$date <= ed
# aux <- mb_weather_all[fil, ]
# save(aux,file = "OntarioSummer2018Hour.RData")

}
\format{A (\eqn{6574 \times 11}{6574 x 11})-matrix containing wind speed observations.}
\source{
  Haslett, J. and Raftery, A. E. (1989),
  Space-time modelling with long-memory dependence: assessing
  Ireland's wind-power resource (with discussion),
  \emph{Applied Statistics}, 38, 1--50.
}

\keyword{datasets}
