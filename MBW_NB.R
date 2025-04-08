#Danielle's code to look at the effect of year on detection. 

#Libraries
library(naturecounts)
library(tidyr)

rm(list=ls())
par(mfrow=c(1,1),mar=c(5,5,3,2)+0.1,cex.axis=1.5,cex.lab=1.5)

###Full MBW dataset from NatureCounts###
#replace with your user name and request_id
#only need to do once per season
MBW.NC <- nc_data_dl(collection="MBW", fields_set = "extended", username = "dethier",
                     info = "analysis for HELP data internal")

#write to data folder
write.csv(MBW.NC, "Data/MBW_2024.csv")

ddf<-read.csv("Data/MBW_2024.csv") #pull in data from file for working with

#retain only the columns that will be useful for the analysis
ddf<-ddf %>% select(SurveyAreaIdentifier, RouteIdentifier, species_id, CommonName, survey_year, survey_month, survey_day, TimeObservationsStarted, TimeObservationsEnded, ObservationCount, ObservationDescriptor, ObservationCount2, ObservationDescriptor2, ObservationCount3, ObservationDescriptor3, ObservationCount4, ObservationDescriptor4, 
                    ObservationCount5, ObservationDescriptor5, ObservationCount6, ObservationDescriptor6, ObservationCount7, ObservationDescriptor7, ObservationCount8, ObservationDescriptor8,  ObservationCount9, ObservationDescriptor9, EffortMeasurement1, EffortUnits1, EffortMeasurement2, EffortUnits2, EffortMeasurement3, EffortUnits3, EffortMeasurement4, EffortUnits4, CollectorNumber, DecimalLatitude, DecimalLongitude)

#create doy field
ddf<-ddf %>% format_dates()

#Need to summarize the data (0-5, 5-10, 10-15, 15-20)
ddf<-ddf %>% rowwise() %>% 
  mutate(Obs1 = sum(ObservationCount2+ObservationCount3, na.rm=TRUE), Obs2= sum(ObservationCount4+ObservationCount5, na.rm=TRUE), Obs3 = sum(ObservationCount6+ObservationCount7, na.rm=TRUE), Obs4=sum(ObservationCount8+ObservationCount9, na.rm=TRUE))

#Create Events Matrix which includes the colvarites of interest
events<-ddf %>% select(SurveyAreaIdentifier, survey_year, survey_month, survey_day, doy, TimeObservationsStarted, EffortMeasurement1, EffortMeasurement3, EffortMeasurement4, DecimalLatitude, DecimalLongitude, CollectorNumber) %>% distinct()

ddf<-ddf %>% filter(!is.na(species_id))
sp.list<-unique(ddf$CommonName)

long_df <- ddf %>%
  pivot_longer(
    cols = c(TimeObservationsStarted),  # Replace with your actual column names
    names_to = "TimeStart",             # New column for original column names
    values_to = "Time"                # New column for corresponding values
  )

#create sp loop to get your Abundance Martix

for(m in 1:length(sp.list)) {
  
  # m<-1 #for testing
  
  sp.ddf<-ddf %>% filter(CommonName==sp.list[m])
  sp.ddf<-sp.ddf %>% select(SurveyAreaIdentifier, survey_year, survey_month, survey_day, TimeObservationsStarted, doy, Obs1, Obs2, Obs3, Obs4)
  sp.ddf<-left_join(events, sp.ddf, by=c("SurveyAreaIdentifier", "survey_year", "survey_month", "survey_day", "doy", "TimeObservationsStarted"))
  sp.ddf<-sp.ddf %>% mutate(Obs1=ifelse(is.na(Obs1), 0, Obs1), Obs2=ifelse(is.na(Obs2), 0, Obs2), Obs3=ifelse(is.na(Obs3), 0, Obs3), Obs4=ifelse(is.na(Obs4), 0, Obs4))
  
  #create species matrix
  sp.matrix<-sp.ddf %>% select(Obs1, Obs2, Obs3, Obs4)
  sp.matrix<-as.matrix(sp.matrix)
  
  write.csv(sp.matrix, paste(sp.list[m], "_DetectMBW.csv", sep=""), row.names = FALSE)
  
  #create covarite matrix
  cov.matrix<-sp.ddf %>% select(survey_year)
  # Copy the column 3 times
  cov.matrix <- cbind(cov.matrix, replicate(3, cov.matrix$survey_year))
  colnames(cov.matrix)[1:4] <- paste0("SampleYear", 1:4)
  cov.matrix<-as.matrix(cov.matrix)
  
  write.csv(cov.matrix, paste(sp.list[m], "_CovMBW.csv", sep=""), row.names = FALSE)
  
} # end sp.loop


##Now we want to assess the effect of year on detection probability

for(m in 1:length(sp.list)) {
  
y<-read.csv(paste(sp.list[m], "_DetectMBW.csv", sep=""))
ObsCovs<-read.csv(paste(sp.list[m], "_CovMBW.csv", sep=""))
#standardize year
ObsCovs <- as.data.frame(lapply(ObsCovs, as.numeric))
ObsCovs <- scale(ObsCovs) 
ObsCovs<-as.matrix(ObsCovs)

# Create unmarkedFrame without site covariates
umf <- unmarkedFramePCount(
  y = y,          # Detection matrix (sites x surveys)
 obsCovs = list(year = ObsCovs) )  # Year covariate matrix (same dim as y)


# Fit model with year affecting detection probability
fm <- pcount(~ year             # Detection formula (year effect)
             ~ 1,              # Abundance formula (intercept only)
             data = umf,
             K = 100)          # Upper integration limit

# Check year's influence on detection
summary(fm)

# Get yearly detection probabilities
det_coef <- coef(fm, type = "det")
std_years <- scale(2012:2024)  # Replace with your actual years
# Predict across standardized years
pred_p$estimate <- plogis(det_coef["p(Int)"] + det_coef["p(year)"] * std_years)

plot<-ggplot(pred_p, aes(x=year.1, y=estimate)) +
  geom_point(size=3) +
  ylim(0,1) +
  labs(title= "Detection Probability by Year", y="p̂")

ggsave(paste(sp.list[m], "DetectionPropbyYear.pdf", sep=""), plot = plot, width = 8, height = 6)

}
