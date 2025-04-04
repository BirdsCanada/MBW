## Mountain Birdwatch Data Import #############################################


# Set up work space + load appropriate packages etc

#setwd(choose.dir())
setwd("C:/Users/akouw/OneDrive/Desktop/R/MBW2024") 

library(dplyr)
library(unmarked)
library(ggplot2)
library(lubridate)
library(MuMIn)
library(naturecounts)

#collection<-meta_collections()

###Full MBW dataset from NatureCounts###
nc_requests(username = "amyleek")
MBW.NC <- nc_data_dl(request_id = 252290, fields_set = "extended", username = "amyleek",
                      info = "analysis for HELP data internal")


#####################
#Zero-filling Dataset
#####################

#########################
#Using zero-fill function
#########################

sp.list <- search_species()
sp.list <- sp.list %>% select(species_id,english_name)

MBW.NC <- MBW.NC %>% left_join(sp.list, by = "species_id")
MBW.NC$ObservationCount <- as.numeric(MBW.NC$ObservationCount)

count(MBW.NC, AllSpeciesReported)

MBW.NC_zero_fill <- format_zero_fill(MBW.NC, 
                                        by = "SamplingEventIdentifier",
                                        extra_event = c("latitude", "longitude"),
                                        extra_species = "english_name")


MBW.NC_presence <- MBW.NC %>%
  select(AllSpeciesReported, SamplingEventIdentifier, species_id, ObservationCount, latitude, longitude, english_name) %>%
  mutate(presence = if_else(as.numeric(ObservationCount) > 0, TRUE, FALSE))


MBW.NC_pres_filled <- format_zero_fill(MBW.NC_presence, 
                                          fill = "presence",
                                          by = "SamplingEventIdentifier",
                                          extra_event = c("latitude", "longitude"), 
                                          extra_species = "english_name")

####################
#Using events matrix
####################
nc_requests(username = "amyleek")
MBW.NC <- nc_data_dl(request_id = 252290, fields_set = "extended", username = "amyleek",
                     info = "analysis for HELP data internal")

events_matrix <- MBW.NC %>%
  select(SurveyAreaIdentifier,
         survey_year,
         survey_month,
         survey_day,
         latitude, 
         longitude) %>%
  distinct()


search_species("Bicknell's Thrush")

bith_events <- MBW.NC %>%
  filter(species_id == "15570")

bith_events <- left_join(events_matrix, bith_events,
                             by = c("SurveyAreaIdentifier",
                                    "survey_year",
                                    "survey_month",
                                    "survey_day",
                                    "latitude",
                                    "longitude"))


bith_zero_fill <- bith_events %>% 
  mutate(ObservationCount = replace(ObservationCount, is.na(ObservationCount), 0)) %>% select(SurveyAreaIdentifier, survey_day, survey_month, survey_year, latitude, longitude, ObservationCount) 


bith_zero_fill$ObservationCount <- as.numeric(bith_zero_fill$ObservationCount)

hist(bith_zero_fill$ObservationCount, breaks = 10)


#########################
##Loop for target species
#########################

nc_requests(username = "amyleek")
MBW.NC <- nc_data_dl(request_id = 252290, fields_set = "extended", username = "amyleek",
                     info = "analysis for HELP data internal")

table(MBW.NC$species_id)
table(MBW.NC$CommonName)

sp_ids <- unique(MBW.NC$species_id)
sp_ids


zero_filled_data <- list()

# Loop through each sp_ids
for (i in sp_ids) {
  # Filter the BBS data for each species in a given loop
  all_species_data <- MBW.NC %>% 
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


####################
####################





######################
####MAPPING####
######################

library(naturecounts)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(ggspatial)
library(dplyr)
library(tidyr)
library(mapview)
library(prettymapr)
library(tidyverse)


#All years

dat<-MBW.NC %>% select(RouteIdentifier, SurveyAreaIdentifier, SamplingEventIdentifier, ProtocolCode, DecimalLatitude, DecimalLongitude, YearCollected, MonthCollected, DayCollected, CommonName, SpeciesCode, Collector, TimeCollected, TimeObservationsStarted, TimeObservationsEnded, EffortMeasurement1, EffortUnits1, EffortMeasurement3, EffortUnits3, EffortMeasurement4, EffortUnits4, ObservationCount, ObservationDescriptor, ObservationCount2, ObservationDescriptor2, ObservationCount3, ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, ObservationCount5, ObservationDescriptor5, ObservationCount6, ObservationDescriptor6)

# JULIAN DATE for time sequence analysis
dat$date<-paste(dat$DayCollected,dat$MonthCollected, dat$YearCollected,sep="/")
dat$date<-as.Date(dat$date, format='%d/%m/%Y')
dat$doy<-as.numeric(format(dat$date, format="%j"))

##Remove NAs and get numerics and factors right
dat <- drop_na(dat, "ObservationCount")
dat$ObservationCount = as.numeric(dat$ObservationCount)
dat$YearCollected <- as.factor(dat$YearCollected)


#Plot raw counts
ggplot(data = dat) + 
  geom_point(aes(x = doy, y = ObservationCount))


#Summarize
#detach(package:plyr)

se <- function(x) sd(x) / sqrt(length(x))

#By year
dat_year <- dat %>% 
  group_by(YearCollected) %>% 
  summarise(MeanObs = mean(ObservationCount), 
            SEObs = se(ObservationCount)) %>%   
  mutate(yrmax = MeanObs + SEObs, yrmin = MeanObs - SEObs)

#Plot
ggplot(data = dat_year) +  
  geom_pointrange(aes(x = YearCollected, y = MeanObs, ymin = yrmin, ymax = yrmax))


#By doy

dat_doy_yr <- dat%>% 
  group_by(YearCollected, doy) %>% 
  summarise(MeanObs = mean(ObservationCount), 
            SEObs = se(ObservationCount)) %>%   
  mutate(yrmax = MeanObs + SEObs, yrmin = MeanObs - SEObs)

#Plot
ggplot(data = dat_doy_yr) +  
  geom_pointrange(aes(x = doy, y = MeanObs, ymin = yrmin, ymax = yrmax, colour = YearCollected))


#Make spatial
dat <- drop_na(dat, "DecimalLongitude", "DecimalLatitude")
dat_sf <- st_as_sf(dat, 
                      coords = c("DecimalLongitude", "DecimalLatitude"), crs = 4326)

ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0) +
  geom_sf(data = dat_sf) +
  labs(caption = "Copyright OpenStreetMap contributors")

#interactive maps
mapview(dat_sf, zcol = "YearCollected", at = seq(2012, 2024, by = 1),
        map.types = "Esri.WorldImagery")


##Making Shapefiles

# Write the shapefile
st_write(dat_sf, "MBW_shapefile.shp")

# Read the shapefile
check_shapefile <- st_read("MBW_shapefile.shp")

# View the first few rows
head(check_shapefile)




##########################################
#### USING JUST BITH DATA ###
##########################################

#nc_requests(username = "amyleek")
#MBW.NC <- nc_data_dl(request_id = 252290, fields_set = "extended", username = "amyleek",
#                     info = "analysis for HELP data internal")

#dat<-MBW.NC %>% select(RouteIdentifier, SurveyAreaIdentifier, SamplingEventIdentifier, ProtocolCode, DecimalLatitude, DecimalLongitude, YearCollected, MonthCollected, DayCollected, CommonName, SpeciesCode, Collector, TimeCollected, TimeObservationsStarted, TimeObservationsEnded, EffortMeasurement1, EffortUnits1, EffortMeasurement3, EffortUnits3, EffortMeasurement4, EffortUnits4, ObservationCount, ObservationDescriptor, ObservationCount2, ObservationDescriptor2, ObservationCount3, ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, ObservationCount5, ObservationDescriptor5, ObservationCount6, ObservationDescriptor6)

dat.BITH <-filter(dat, SpeciesCode == "15570")
#dat.24 <-filter(dat, YearCollected == "2024")
#dat.BITH24 <-filter(dat.BITH, YearCollected == "2024")



ggplot(data = dat.BITH) + 
  geom_point(aes(x = doy, y = ObservationCount))


#Summarize
#detach(package:plyr)

se <- function(x) sd(x) / sqrt(length(x))

#By year
dat.BITH_year <- dat.BITH %>% 
  group_by(YearCollected) %>% 
  summarise(MeanObs = mean(ObservationCount), 
            SEObs = se(ObservationCount)) %>%   
  mutate(yrmax = MeanObs + SEObs, yrmin = MeanObs - SEObs)

#Plot
ggplot(data = dat.BITH_year) +  
  geom_pointrange(aes(x = YearCollected, y = MeanObs, ymin = yrmin, ymax = yrmax))


#By doy

dat.BITH_doy_yr <- dat.BITH %>% 
  group_by(YearCollected, doy) %>% 
  summarise(MeanObs = mean(ObservationCount), 
            SEObs = se(ObservationCount)) %>%   
  mutate(yrmax = MeanObs + SEObs, yrmin = MeanObs - SEObs)

#Plot
ggplot(data = dat.BITH_doy_yr) +  
  geom_pointrange(aes(x = doy, y = MeanObs, ymin = yrmin, ymax = yrmax, colour = YearCollected))


###Mapping just BITH

#Make spatial
dat.BITH <- drop_na(dat.BITH, "DecimalLongitude", "DecimalLatitude")
dat.BITH_sf <- st_as_sf(dat.BITH, 
                        coords = c("DecimalLongitude", "DecimalLatitude"), crs = 4326)

ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0) +
  geom_sf(data = dat.BITH_sf) +
  labs(caption = "Copyright OpenStreetMap contributors")

#interactive maps
mapview(dat.BITH_sf, zcol = "YearCollected", at = seq(2012, 2024, by = 1),
        map.types = "Esri.WorldImagery")


##Making Shapefiles

# Write the shapefile
st_write(dat.BITH_sf, "BITH_shapefile2.shp")

# Read the shapefile
check_shapefile <- st_read("BITH_shapefile2.shp")

# View the first few rows
head(check_shapefile)



 
