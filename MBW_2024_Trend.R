#Libraries
library(naturecounts)
library(tidyr)
library(dplyr)
library(unmarked)
library(ggplot2)
library(lubridate)
library(MuMIn)
library(DHARMa)


rm(list=ls())
par(mfrow=c(1,1),mar=c(5,5,3,2)+0.1,cex.axis=1.5,cex.lab=1.5)

###Full MBW dataset from NatureCounts###
#replace with your user name and request_id
#only need to do once per season
nc_requests(username = "amyleek")
MBW.NC <- nc_data_dl(request_id = 252290, fields_set = "extended", username = "amyleek",
                     info = "analysis for HELP data internal")

#write to data folder
write.csv(MBW.NC, "Data/MBW_2024.csv")

ddf<-read.csv("Data/MBW_2024.csv") #pull in data from file for working with

#remove NS data from 2016 and 2017
ddf <- filter(ddf, subnational2_code != "CA.NS.IN", subnational2_code != "CA.NS.VI")


#########################
##Loop for target species
#########################

events_matrix <- ddf %>%
  select(SurveyAreaIdentifier,
         survey_year,
         survey_month,
         survey_day,
         latitude, 
         longitude) %>%
  distinct()

table(ddf$species_id)
table(ddf$CommonName)

sp_ids <- unique(ddf$species_id)
sp_ids


zero_filled_data <- list()

# Loop through each sp_ids
for (i in sp_ids) {
  # Filter the BBS data for each species in a given loop
  all_species_data <- ddf %>% 
    filter(species_id == i)
  
  # Join the NatureCounts data with the events matrix
  all_species_events <- left_join(events_matrix, all_species_data,
                                  by = c("SurveyAreaIdentifier", 
                                         "survey_year", 
                                         "survey_month",
                                         "survey_day", 
                                         "latitude", 
                                         "longitude"))
  
  # Zero-fill the `ObservationCount` column for each species and across all events
  all_species_events <- all_species_events %>%
    mutate(ObservationCount = replace(ObservationCount, is.na(ObservationCount), 0))
  
  # Add the result to the list with species_id as the key
  zero_filled_data[[as.character(i)]] <- all_species_events
}

zero_filled_dataframe <- bind_rows(zero_filled_data, .id = "species_id")


ddf.zf <- zero_filled_dataframe




#retain only the columns that will be useful for the analysis
ddf.zf<-ddf.zf %>% select(SurveyAreaIdentifier, RouteIdentifier, ProtocolCode, species_id, CommonName, subnational2_code, survey_year, survey_month, survey_day, TimeObservationsStarted, TimeObservationsEnded, ObservationCount, ObservationDescriptor, ObservationCount2, ObservationDescriptor2, ObservationCount3, ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, 
                    ObservationCount5, ObservationDescriptor5, ObservationCount6, ObservationDescriptor6, ObservationCount7, ObservationDescriptor7, ObservationCount8, ObservationDescriptor8,  ObservationCount9, ObservationDescriptor9, EffortMeasurement1, EffortUnits1, EffortMeasurement2, EffortUnits2, EffortMeasurement3, EffortUnits3, EffortMeasurement4, EffortUnits4, CollectorNumber, DecimalLatitude, DecimalLongitude, AllSpeciesReported)
#remove NS data from 2016 and 2017
#ddf <- filter(ddf, subnational2_code != "CA.NS.IN", subnational2_code != "CA.NS.VI")

#create doy field
ddf.zf<-ddf.zf %>% format_dates()


####NEW ANALYSIS###

install.packages(c("performance", "see"))
library(performance)
library(see)
library(lme4)

#make year continous
ddf.zf$survey_year <- as.numeric(ddf.zf$survey_year)

#Calculate survey start column. Time since sunrise (~5:30AM sunrise Popple Depot, NB in June, but used 4:15AM to avoid negative numbers.)

ddf.zf$TimeSinceSunrise <- ddf.zf$TimeObservationsStarted - 4.25
#ddf.zf$TimeBeforeSunset <- 9.5- ddf.zf$TimeObservationsStarted

ddf.zf.bith <- subset(ddf.zf, CommonName == "Bicknell's Thrush")
#ddf.zf.wtsp <- subset(ddf.zf, CommonName == "White-throated Sparrow")


#Poisson GLM

#GLM1 <- glm(ObservationCount ~ doy + survey_year + TimeSinceSunrise, family = poisson, data = ddf.zf.bith)

#GLM2 <- glmer(ObservationCount ~ doy + survey_year + TimeSinceSunrise + (1 | RouteIdentifier),
#             family = poisson, data = ddf.zf.bith)

library(glmmTMB)
GLM3 <- glmmTMB(ObservationCount ~ survey_year + doy + TimeSinceSunrise + (1 | RouteIdentifier),
              family = poisson, data = ddf.zf.bith)

# Generate diagnostic plots
check_model(GLM3)
summary(GLM3)



plot(ObservationCount ~ TimeSinceSunrise, data = ddf.zf.bith)
plot(ObservationCount ~ doy, data = ddf.zf.bith)
plot(ObservationCount ~ survey_year, data = ddf.zf.bith)

