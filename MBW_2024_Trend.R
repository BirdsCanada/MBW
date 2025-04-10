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

#get rid of survey with weird start time
ddf <- filter(ddf, TimeObservationsStarted != 13.3333)


#retain only the columns that will be useful for the analysis
ddf<-ddf %>% select(SurveyAreaIdentifier, RouteIdentifier, ProtocolCode, species_id, CommonName, subnational2_code, survey_year, survey_month, survey_day, TimeObservationsStarted, TimeObservationsEnded, ObservationCount, ObservationDescriptor, ObservationCount2, ObservationDescriptor2, ObservationCount3, ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, 
                    ObservationCount5, ObservationDescriptor5, ObservationCount6, ObservationDescriptor6, ObservationCount7, ObservationDescriptor7, ObservationCount8, ObservationDescriptor8,  ObservationCount9, ObservationDescriptor9, EffortMeasurement1, EffortUnits1, EffortMeasurement2, EffortUnits2, EffortMeasurement3, EffortUnits3, EffortMeasurement4, EffortUnits4, CollectorNumber, DecimalLatitude, DecimalLongitude, AllSpeciesReported)

#create doy field
ddf<-ddf %>% format_dates()


####NEW ANALYSIS###

#install.packages(c("performance", "see"))
library(performance)
library(see)
library(lme4)

#make year continous
ddf$survey_year <- as.numeric(ddf$survey_year)

#Calculate survey start column. Time since sunrise (~5:30AM sunrise Popple Depot, NB in June, but used 4:15AM to avoid negative numbers.)

#ddf$TimeSinceSunrise <- ddf$TimeObservationsStarted - 4.25
#ddf$TimeBeforeSunset <- ddf$TimeObservationsStarted - 20

library(suncalc)
library(ARUtools)


#Poisson GLM

#GLM1 <- glm(ObservationCount ~ doy + survey_year + TimeSinceSunrise, family = poisson, data = ddf.zf.bith)

GLM2 <- glmer(ObservationCount ~ doy + survey_year + ProtocolCode + (1 | RouteIdentifier),
             family = poisson, data = ddf)


library(glmmTMB)
GLM3 <- glmmTMB(ObservationCount ~ survey_year + doy + ProtocolCode + (1 | RouteIdentifier),
              family = poisson, data = ddf)

# Generate diagnostic plots
check_model(GLM2)
summary(GLM2)




###April 9 from Danielle###
#####Loop to zero fill and do analysis####

results <- data.frame(Group = integer(),
                      Estimate = numeric(),
                      Std.Error = numeric(),
                      z.value = numeric(),
                      p.value = numeric(),
                      stringsAsFactors = FALSE)


sp_ids <- unique(ddf$species_id)

for(m in 1:length(sp_ids)) {
  
   #m<-1 #for testing
  
  sp.ddf<-ddf %>% filter(CommonName==sp_ids[m]) #this will cycle through each species in sp.ids. For testing you can manually set m to 
  sp.ddf<-sp.ddf %>% select(SurveyAreaIdentifier, RouteIdentifier, survey_year, doy, TimeObservationsStarted)
  sp.ddf<-left_join(all_species_events, sp.ddf, by=c("SurveyAreaIdentifier", "survey_year")) #you will zero fill in the loop
    
              
                    ###now you will put your analysis code here
                    ##Then you will put your output table script here

  GLM2 <- glmer(ObservationCount ~ doy + survey_year + ProtocolCode + (1 | RouteIdentifier),
                family = poisson, data = ddf)
  
  # Get summary of the model
  model_summary <- summary(GLM2)
  
  # Extract coefficients and statistics
  estimate <- model_summary$coefficients[2, 1]  # Estimate for predictor
  std_error <- model_summary$coefficients[2, 2]  # Standard error
  z_value <- model_summary$coefficients[2, 3]    # z-value
  p_value <- model_summary$coefficients[2, 4]    # p-value
  
  # Append results to the results data frame
  results <- rbind(results, data.frame(
    Group = sp_ids[m],
    Estimate = estimate,
    Std.Error = std_error,
    z.value = z_value,
    p.value = p_value
  ))
}

print(results)


# This closes the loop

check_model(GLM2)
