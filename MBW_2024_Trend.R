##Libraries

library(naturecounts)
library(tidyr)
library(unmarked)
library(MuMIn)
library(DHARMa)
library(performance)
library(see)
library(lme4)
library(sf)
library(ggspatial)
library(suncalc)
library(ARUtools)
library(glmmTMB)
library(lme4)    

##Load Data

#Full MBW dataset from NatureCounts###
#replace with your user name and request_id
#only need to do once per season
nc_requests(username = "amyleek")
MBW.NC <- nc_data_dl(request_id = 252290, fields_set = "extended", username = "amyleek",
                     info = "analysis for HELP data internal")
#write to data folder
write.csv(MBW.NC, "Data/MBW_2024.csv")

#read from data folder
ddf<-read.csv("Data/MBW_2024.csv") #pull in data from file for working with

##Data Cleaning

#remove NS data from 2016 and 2017
ddf <- filter(ddf, subnational2_code != "CA.NS.IN", subnational2_code != "CA.NS.VI")

#get rid of survey with weird start time
ddf <- filter(ddf, TimeObservationsStarted != 13.3333)

#retain only the columns that will be useful for the analysis
ddf<-ddf %>% select(SurveyAreaIdentifier, RouteIdentifier, ProtocolCode, species_id, CommonName, subnational2_code, survey_year, survey_month, survey_day, TimeObservationsStarted, TimeObservationsEnded, ObservationCount, ObservationDescriptor, ObservationCount2, ObservationDescriptor2, ObservationCount3, ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, 
                    ObservationCount5, ObservationDescriptor5, ObservationCount6, ObservationDescriptor6, ObservationCount7, ObservationDescriptor7, ObservationCount8, ObservationDescriptor8,  ObservationCount9, ObservationDescriptor9, EffortMeasurement1, EffortUnits1, EffortMeasurement2, EffortUnits2, EffortMeasurement3, EffortUnits3, EffortMeasurement4, EffortUnits4, CollectorNumber, DecimalLatitude, DecimalLongitude, AllSpeciesReported)

#create doy field
ddf<-ddf %>% format_dates()

##Data Visualization
ddf_sf <- st_as_sf(ddf, coords = c("DecimalLatitude", "DecimalLongitude"), crs = 4326)

  plot<-ggplot(data = ddf_sf) +
    # Select a basemap
    annotation_map_tile(type = "cartolight", zoom = NULL, progress = "none") +
    # Plot the points, color-coded by survey_year
    geom_sf(aes(color = as.factor(survey_year)), size = 1) +
    # Add a theme with a minimal design and change the font styles, to your preference
    theme_minimal() +
       theme(legend.position = "bottom") +
    # To make the points in the legend larger without affecting map points
    guides(color = guide_legend(override.aes = list(size = 3))) +
    # Define the title and axis names
    labs(title = "MBW Survey Stops",
         x = "Longitude",
         y = "Latitude")
###Need to check your latitude on some of your points. 
  
##Data Summary
  
#total count of species per year
sum_sp<-ddf %>% group_by(CommonName, survey_year) %>% summarise(CountTot=sum(ObservationCount, na.rm =TRUE)) %>% filter(!is.na(CommonName))
sum_sp<-pivot_wider(
  data = sum_sp,
  names_from = survey_year,    # Column to use for new column names
  values_from = CountTot      # Column containing values to fill
)

#Events
#Create Events Matrix which includes the colvarites of interest
all_species_events<-ddf %>% select(SurveyAreaIdentifier, survey_year, survey_month, survey_day, doy, TimeObservationsStarted, DecimalLatitude, DecimalLongitude, ProtocolCode, RouteIdentifier) %>% distinct()

####NEW ANALYSIS###

#make year continuous
ddf$survey_year <- as.numeric(ddf$survey_year)

#Calculate survey start column. Time since sunrise (~5:30AM sunrise Popple Depot, NB in June, but used 4:15AM to avoid negative numbers.)

#ddf$TimeSinceSunrise <- ddf$TimeObservationsStarted - 4.25
#ddf$TimeBeforeSunset <- ddf$TimeObservationsStarted - 20


#Poisson GLM

#GLM1 <- glm(ObservationCount ~ doy + survey_year + TimeSinceSunrise, family = poisson, data = ddf.zf.bith)

GLM2 <- glmer(ObservationCount ~ doy + survey_year + ProtocolCode + (1 | RouteIdentifier),
             family = poisson, data = ddf)

#Danielle new suggested model to try
GLM2 <- glmer(ObservationCount ~  survey_year + ProtocolCode + (1 | SurveyAreaIndentifier),
              family = poisson, data = ddf)


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


sp_ids<-unique(ddf$CommonName)

for(m in 1:length(sp_ids)) {
  
   #m<-1 #for testing
  
  sp.ddf<-ddf %>% filter(CommonName==sp_ids[m]) #this will cycle through each species in sp.ids. For testing you can manually set m to 
  sp.ddf<-sp.ddf %>% select(CommonName, SurveyAreaIdentifier, survey_year, doy, TimeObservationsStarted, ObservationCount)
  sp.ddf<-left_join(all_species_events, sp.ddf, by=c("SurveyAreaIdentifier", "survey_year", "doy", "TimeObservationsStarted")) #you will zero fill in the loop
  #Add the 0 to observation count
  sp.ddf <- sp.ddf %>%
    mutate(ObservationCount = replace(ObservationCount, is.na(ObservationCount), 0))  
              
  hist(sp.ddf$ObservationCount)
  #Prepare variable<-
  sp.ddf$SurveyAreaIdentifier<-as.numeric(factor(paste(sp.ddf$SurveyAreaIdentifier)))
  sp.ddf$ProtocolCode<-as.numeric(factor(paste(sp.ddf$ProtocolCode)))
  sp.ddf$scaleyear<-scale(sp.ddf$survey_year, center = TRUE, scale = TRUE)
  
  GLM1<- glmer(ObservationCount ~ scaleyear + (1 | SurveyAreaIdentifier), data = sp.ddf, family = poisson)
  GLM2 <- glmer.nb(ObservationCount ~ scaleyear + (1 | SurveyAreaIdentifier), data = sp.ddf)
  
  # Get summary of the model
  model_summary1 <-summary(GLM1)
  model_summary2 <- summary(GLM2)
  
  ###I have not done anything past getting the models to run###
  
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
