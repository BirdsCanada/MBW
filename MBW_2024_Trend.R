##Libraries

library(naturecounts)
library(tidyverse)
library(DHARMa)
#library(performance)
#library(see)
library(lme4)
library(sf)
library(ggspatial)
#library(suncalc)
library(ARUtools)
library(glmmTMB)
library(raster)
library(ggmap)
library(patchwork)
library(performance)
library(DHARMa)
library(scales)
   
output_dir <- "Output"
##Load Data

#Full MBW dataset from NatureCounts###
#replace with your user name and request_id
#only need to do once per season
nc_requests(username = "amyleek")
MBW.NC <- nc_data_dl(request_id = 252290, fields_set = "extended", username = "amyleek",
                     info = "analysis for HELP data internal")
#ddf <-MBW.NC
#write to data folder
write.csv(MBW.NC, "Data/MBW_2025_Nov4.csv")

#read from data folder
#ddf<-read.csv("Data/MBW_2025_Nov4_withTotalBITH_NB.csv") #pull in data from file for working with
ddf<-read.csv("Data/MBW_2025_Nov4.csv") #pull in data from file for working with

##Data Cleaning

#remove NS data from 2016 and 2017
ddf <- ddf %>% dplyr::filter(subnational2_code != "CA.NS.IN", subnational2_code != "CA.NS.VI")

#remove prior to 2016
ddf <- ddf %>%  dplyr::filter(survey_year>=2016)

#get rid of survey with weird start time
ddf <- ddf %>% dplyr::filter(TimeObservationsStarted != 16.7167)

#get rid of no observations rows
#ddf <- ddf %>% dplyr::filter(NoObservations != "NoObs" | is.na(NoObservations))

#take out routes only run twice or less in 10 years
#ddf <- ddf %>% dplyr::filter(RouteIdentifier != "NBMBW35", RouteIdentifier != "NBMBW52", RouteIdentifier != "NBMBW65", RouteIdentifier != "NBMBW75",  RouteIdentifier != "NBMBW79", RouteIdentifier != "NBMBW80", RouteIdentifier != "NBMBW97", RouteIdentifier != "NBMBW60", RouteIdentifier != "NBMBW68")


#retain only the columns that will be useful for the analysis
ddf<-ddf %>% dplyr::select(SurveyAreaIdentifier, RouteIdentifier, ProtocolCode, species_id, CommonName, subnational2_code, survey_year, survey_month, survey_day, TimeObservationsStarted, TimeObservationsEnded, ObservationCount, ObservationDescriptor, ObservationCount2, ObservationDescriptor2, ObservationCount3, ObservationDescriptor3, ObservationCount4, ObservationDescriptor4,
                    ObservationCount5, ObservationDescriptor5, ObservationCount6, ObservationDescriptor6, ObservationCount7, ObservationDescriptor7, ObservationCount8, ObservationDescriptor8,  ObservationCount9, ObservationDescriptor9, EffortMeasurement1, EffortUnits1, EffortMeasurement2, EffortUnits2, EffortMeasurement3, EffortUnits3, EffortMeasurement4, EffortUnits4, CollectorNumber, DecimalLatitude, DecimalLongitude, AllSpeciesReported)


#ddf_BITH <- subset(ddf, ddf$CommonName == "Bicknell's Thrush")
#ddf_BITH25 <- subset(ddf_BITH, ddf_BITH$survey_year == 2025)
#write.csv(ddf_BITH25, "BITH_2025.csv")
#write.csv(ddf_BITH, "BITH_allyears.csv")
  
#create doy field
ddf<-ddf %>% format_dates()

##Data Visualization
##Added backgroun map and facet wrap per year to better visualize 
ddf_sf <- st_as_sf(ddf, coords = c("DecimalLatitude", "DecimalLongitude"), crs = 4326)

# Create the plot
plot <- ggplot2::ggplot(data = ddf_sf) +
  # Select a basemap
  annotation_map_tile(type = "cartolight", zoom = NULL, progress = "none") +
  # Plot the points, color-coded by survey_year
  geom_sf(aes(color = as.factor(survey_year)), size = 1) +
  # Add the facet wrap to create a panel for each year
  facet_wrap(~ survey_year) +
  # Add a theme with a minimal design
  theme_minimal() +
  theme(legend.position = "bottom") +
  # To make the points in the legend larger without affecting map points
  guides(color = guide_legend(override.aes = list(size = 3))) +
  # Define the title and axis names
  labs(title = "MBW Survey Stops by Year",
       x = "Longitude",
       y = "Latitude",
       color = "Survey Year") # Adds a title to the legend

# Display the plot This is ugly but shows everything looks OK. 
plot

ddf$ObservationCount<-as.integer(ddf$ObservationCount)
#total count of species per year
sum_sp1_stats <- ddf %>%
  filter(!is.na(CommonName)) %>% 
  group_by(CommonName, survey_year) %>%
  summarise(
    # The total count you already had
    CountTot = sum(ObservationCount, na.rm = TRUE),
    # Calculate Mean, Median, Min, and Max for each group
    MeanCount = mean(ObservationCount, na.rm = TRUE),
    MedianCount = median(ObservationCount, na.rm = TRUE),
    MinCount = min(ObservationCount, na.rm = TRUE),
    MaxCount = max(ObservationCount, na.rm = TRUE),
    # It's good practice to also count the number of observations in each group
    N = n(),
    .groups = 'drop' # This is good practice to prevent issues later
  )

sum_sp1<- sum_sp1_stats %>% dplyr::select(CommonName, survey_year, CountTot)

sum_sp<-pivot_wider(
  data = sum_sp1,
  names_from = survey_year,    # Column to use for new column names
  values_from = CountTot      # Column containing values to fill
)

write.csv(sum_sp, "TotalCountSpeciesPerYear25_4Nov.csv")
ggplot(data = sum_sp1)+ 
  geom_point(aes(x = survey_year, y = CountTot))

ggplot(data = sum_sp1) +  
  geom_pointrange(aes(x = survey_year, y = CountTot, ymin = 0, ymax = 600, colour = CommonName))


#Multiple lines and plots
#sum_sp1$survey_year <- as.numeric(sum_sp1$survey_year)

ggplot(sum_sp1, aes(survey_year, CountTot, colour = CommonName)) +  
  geom_point()+
  geom_smooth(formula = y ~ x, method = "lm")+
  facet_wrap(~CommonName, scales = "free")+ #set scale = free to better see the differences
  scale_x_continuous(breaks = pretty_breaks())+
labs(
  x = "Survey Year",
  y = "Individuals Detected"
) +
  theme(
    legend.position = "none" # The legend is redundant because of the facet titles
  )
ggsave("Total Count per Year by Species_4Nov.pdf", width = 11, height = 8.5, units ="in")


# Create the box and whisker plot
ggplot(ddf, aes(x = as.factor(survey_year), y = ObservationCount, fill = CommonName)) +
  geom_boxplot() +
  # Use facet_wrap to create a separate panel for each species (CommonName)
  facet_wrap(~ CommonName, scales = "free_y") +
  # Add some labels and a theme for readability
  labs(
    title = "Observation Counts by Year and Species",
    x = "Survey Year",
    y = "Observation Count"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none" # The legend is redundant because of the facet titles
  )


#more basic plot
#sum_sp1$CountTot <- as.numeric(sum_sp1$CountTot)
#sum_sp1$CommonName <- as.factor(sum_sp1$CommonName)


#png("Detections by species and year_31Oct.png")
#par(mfrow = c(2,5))
#plot(CountTot ~ survey_year ,data=sum_sp1, type="p",col = CommonName, ylim=c(0,550),lwd=2,ylab="Number of Individuals Detected",las=2,xlab="",main="Total Individuals Detected by Species and Year",bty="l", xaxt = "n")
#axis(1, at=seq(2016,2025,by = 1), las = 2)

#points(CountTot ~ CommonName, pch=20,col= survey_year,data=sum_sp1)
#lines(CountTot~survey_year, type = "b", col = CommonName, data = sum_sp1)
#dev.off()

#make year continuous
is.numeric(ddf$survey_year)
#ddf$survey_year <- as.numeric(ddf$survey_year)

#Route level Effort
#Determine the number of stops per route to see if we need to effort correct
stops<-ddf %>% group_by(RouteIdentifier, survey_year) %>% summarise(nstop = n_distinct(SurveyAreaIdentifier)) #we will include this as an offset in the model

#Add in route start time with slice_min
start<-ddf %>% dplyr::select(RouteIdentifier, survey_year, TimeObservationsStarted) %>% group_by(RouteIdentifier, survey_year) %>% slice_min(TimeObservationsStarted) %>% distinct()

##Add Red Squirrel for Winter Wren
RS<-sum_sp1 %>% filter(CommonName=="North American Red Squirrel") %>% rename( RedSquirrel = CountTot ) %>% dplyr::select(-CommonName)

#Route level Events
#Create Events Matrix which includes the colvarites of interest
#Removed Survey Area Identifier since we will do the analysis at the route level
all_species_events<-NULL
all_species_events<-ddf %>% dplyr::select(survey_year, survey_month, survey_day, doy, ProtocolCode, RouteIdentifier, SurveyAreaIdentifier, CollectorNumber) %>% distinct()
all_species_events<-all_species_events %>% left_join(stops, by=c("RouteIdentifier", "survey_year"))
all_species_events<-all_species_events %>% left_join(start, by=c("RouteIdentifier", "survey_year"))
all_species_events<-all_species_events %>% left_join(RS, by=c("survey_year"))

results<-NULL #clear old
results <- data.frame(Group = integer(),
                      Estimate = numeric(),
                      Std.Error = numeric(),
                      z.value = numeric(),
                      p.value = numeric(),
                      Dispersion.ratio = numeric(), 
                      Dispersion_p.value = numeric(),
                      Percent_Change = numeric(), 
                      stringsAsFactors = FALSE)


ddf<-ddf %>% filter(!is.na(CommonName))
#ddf<-ddf %>% filter(CommonName != "North American Red Squirrel", CommonName != "Black-capped Chickadee")
sp_ids<-unique(ddf$CommonName)


sp_ids <- sp_ids[!sp_ids %in% c("Black-capped Chickadee", "North American Red Squirrel")]


for(m in 1:length(sp_ids)) {
  
  #m<-3 #for testing
  
  sp.ddf<-NULL #clear old dataframe
  sp.ddf<-ddf %>% filter(CommonName==sp_ids[m]) #this will cycle through each species in sp.ids. For testing you can manually set m to 
  sp.ddf<-sp.ddf %>% dplyr::select(CommonName, RouteIdentifier, SurveyAreaIdentifier, survey_year, doy, ObservationCount)
  sp.ddf<-left_join(all_species_events, sp.ddf, by=c("RouteIdentifier", "SurveyAreaIdentifier", "survey_year", "doy")) #you will zero fill in the loop
  #Add the 0 to observation count
  sp.ddf <- sp.ddf %>%
   mutate(ObservationCount = replace(ObservationCount, is.na(ObservationCount), 0))  
  
  #Sum the total count of individual per route as response
  sp.ddf<-sp.ddf %>% group_by(RouteIdentifier, survey_year, ProtocolCode, nstop, doy, CollectorNumber, RedSquirrel) %>% summarise(RouteTotal = sum(ObservationCount, na.rm=TRUE))
  
  # #Make response variable binomial
  #library(dplyr)
  sp.ddf <- sp.ddf %>%
    mutate(RouteTotal_nb = ifelse(RouteTotal != 0, 1, 0))
  
  
  #Prepare variable
  sp.ddf$RouteIdentifierFact<-as.numeric(factor(paste(sp.ddf$RouteIdentifier)))
  sp.ddf$ProtocolCode<-as.numeric(factor(paste(sp.ddf$ProtocolCode)))
  sp.ddf$CollectorNumber<-as.numeric(factor(paste(sp.ddf$CollectorNumber)))
  sp.ddf$scaleyear<-scale(sp.ddf$survey_year, center = TRUE, scale = TRUE)
  sp.ddf$scalesquirrel<-scale(sp.ddf$RedSquirrel, center = TRUE, scale = TRUE)
  
 
  #remove routes where a species was never detected, or if the route was sampled in < 2 years
  route_remove<-sp.ddf %>% group_by(RouteIdentifier) %>% summarise(keep_route = sum(RouteTotal, na.rm = TRUE), samples = length(RouteTotal)) %>% filter(keep_route>0 & samples >2) %>% dplyr::select(-keep_route, -samples)
  sp.ddf<-left_join(route_remove, sp.ddf, by="RouteIdentifier")

  #create a species specific summary 
  sp.sum<-sp.ddf %>% group_by(survey_year, RouteIdentifier) %>% summarise(n= sum(RouteTotal))
  
  sp_wide <- sp.sum %>%
    pivot_wider(
      names_from = survey_year,     # column headers
      values_from = n              # cell values
    )
  
  write.table(
    sp_wide,
    file = file.path(output_dir, paste0(sp_ids[m], "SpeciesByYearSummary.csv")),
    sep = ",",
    row.names = FALSE
  )

  
  #Sum the count on a given route
  png(filename = file.path(output_dir, paste0(sp_ids[m], "RouteTotal_histogram.png")))
  hist(sp.ddf$RouteTotal)
  dev.off()
  

  
  if(sp_ids[m] %in% c("Hermit Thrush", "Fox Sparrow")){
    
    GLM<- glmmTMB(RouteTotal ~ scaleyear +  ProtocolCode + (1 | RouteIdentifierFact) + offset(log(nstop)), data = sp.ddf, family = nbinom1())
    
  }

 
  if(sp_ids[m] %in% c("White-throated Sparrow")){
    
    GLM<- glmmTMB(RouteTotal ~ scaleyear +  ProtocolCode + (1 | RouteIdentifierFact) + offset(log(nstop)), data = sp.ddf, family = genpois(link = "log"))               
    
  }
  
  
  if(sp_ids[m] %in% c("Blackpoll Warbler")){
    
    GLM<- glmmTMB(RouteTotal ~ scaleyear +  ProtocolCode + (1 | RouteIdentifierFact) + offset(log(nstop)), data = sp.ddf, family = genpois(link = "log"))               
    
  }
  
  if(sp_ids[m] == "Winter Wren"){
    
    GLM<- glmmTMB(RouteTotal ~ scaleyear +  ProtocolCode + scalesquirrel + (1 | RouteIdentifierFact) + offset(log(nstop)), data = sp.ddf, family = genpois(link = "log"))               
    
  }


 
  if(sp_ids[m] %in% c("Bicknell's Thrush", "Yellow-bellied Flycatcher", "Boreal Chickadee")){
    
    GLM <- glmmTMB(RouteTotal ~ scaleyear + ProtocolCode + (1 | RouteIdentifierFact) + offset(log(nstop)), data = sp.ddf,
                   
                   ziformula = ~ 1,  # Intercept-only zero-inflation
                   
                   family = genpois(link = "log"))
                   
  }  
    
    if(sp_ids[m] %in% c("Swainson's Thrush")){
      
      
      GLM <- glmmTMB(RouteTotal ~ scaleyear + ProtocolCode + (1 | RouteIdentifierFact) + offset(log(nstop)), data = sp.ddf,
                     
                     
                     family = gaussian(link = "identity"))  # Use normal distribution
                     
                     
  }

  
  model_summary <- summary(GLM)
 
  # Simulate residuals

  simulationOutput <- simulateResiduals(fittedModel = GLM)


  #QQ plot on the left should follow the line
  #Residual on the right should be random 
  
  # Plot the main diagnostic plot
  pdf_path <- file.path(output_dir, paste(sp_ids[m], "_DHARMa_diagnostics.pdf"))
  pdf(file = pdf_path, width = 7, height = 7)
  plot(simulationOutput)
  plotResiduals(simulationOutput, form = sp.ddf$scaleyear)
  plotResiduals(simulationOutput, form = sp.ddf$ProtocolCode)
  dev.off()

  
  #Test for Over or Under Disperson
  dispersion_result <- testDispersion(simulationOutput)
  dispersion_ratio <- unname(dispersion_result$statistic)
  dispersion_p_value <- dispersion_result$p.value
  #significant result would suggest that the nbinom2 is not appropriate. 
  
  # Access the table of coefficients for the conditional model
  coeffs_table <- model_summary$coefficients$cond
  
    # Extract values for the first predictor (row 2)
  # Column 1: Estimate
  # Column 2: Standard Error
  # Column 3: z-value
  # Column 4: p-value
  estimate <- coeffs_table[2, 1]
  std_error <- coeffs_table[2, 2]
  z_value <- coeffs_table[2, 3]
  p_value <- coeffs_table[2, 4]
  
  # Get standard deviation of original year if scaleyear is scaled
  sd_year <- sd(sp.ddf$survey_year)
  
  # Convert to per-year effect
  beta_year <- estimate / sd_year
  
  # Population change per year (multiplicative factor)
  annual_change_factor <- exp(beta_year)
  
  # Percent change per year
  annual_percent_change <- (annual_change_factor - 1) * 100

  # Append results to the results data frame
  results <- rbind(results, data.frame(
    Group = sp_ids[m],
    Estimate = estimate,
    Std.Error = std_error,
    z.value = z_value,
    p.value = p_value,
    Dispersion.ratio = dispersion_ratio, 
    Dispersion_p.value = dispersion_p_value,
    Annual_Per_Change = annual_percent_change
    
  ))

  
}

file_path2 <- file.path(output_dir, paste0("ModelResults.csv"))
write.table(results, file = file_path2, row.names = FALSE, sep = ",")



# This closes the loop

##################################

####Red Squirrel vs Winter Wren###

##################################

# Aggregate data by year
total_counts_RS <- ddf %>%
  group_by(survey_year) %>%
  summarise(
    total_RS = sum(ObservationCount[CommonName == "North American Red Squirrel"], na.rm = TRUE),
    total_WW = sum(ObservationCount[CommonName == "Winter Wren"], na.rm = TRUE),
    total_BT = sum(ObservationCount[CommonName == "Bicknell's Thrush"], na.rm = TRUE),
    total_ST = sum(ObservationCount[CommonName == "Swainson's Thrush"], na.rm = TRUE),
    total_BP = sum(ObservationCount[CommonName == "Blackpoll Warbler"], na.rm = TRUE)
  )

# Calculate the correlation coefficient
correlation_coefficient <- cor(total_counts_RS$total_RS, total_counts_RS$total_WW)
print(correlation_coefficient)

# Test for significance
correlation_test <- cor.test(total_counts_RS$total_RS, total_counts_RS$total_WW)
print(correlation_test)

# Visualizing the results
colors <- rainbow(length(unique(total_counts_RS$survey_year)))

plot(total_counts_RS$total_RS, total_counts_RS$total_WW,
     main="Comparison of Species Totals Each Year",
     xlab="Total Red Squirrel",
     ylab="Total Winter Wren",
    pch = 19,  # Solid points
    col = colors[as.numeric(factor(total_counts_RS$survey_year))],  # Color by year
abline(lm(total_counts_RS$total_RS ~ total_counts_RS$total_WW), col="red"))

# Add a legend
legend("topright", legend = unique(total_counts_RS$survey_year), 
       col = colors, pch = 19, cex = 0.3, title = "Year")


###Bird Species at year lag compared to Red Squirrel###

# Create a lagged variable for Species A
total_counts_RS <- total_counts_RS %>%
  mutate(total_RS_lag = lag(total_RS, 1)) %>%
  na.omit()  # Remove NA values from lagging

# Check the final dataset
print(total_counts_RS)

# Calculate the correlation coefficient
correlation_coefficient_lag <- cor(total_counts_RS$total_RS_lag, total_counts_RS$total_WW)
print(correlation_coefficient_lag)

# Test for significance
correlation_test_lag <- cor.test(total_counts_RS$total_RS_lag, total_counts_RS$total_WW)
print(correlation_test_lag)

# Visualizing the results
plot(total_counts_RS$total_RS_lag, total_counts_RS$total_WW, 
     main="Lagged Comparison of Species Counts",
     xlab="Winter Wren",
     ylab="Total Red Squirrel (Lagged)")
abline(lm(total_counts_RS$total_RS_lag ~ total_counts_RS$total_WW), col="blue")