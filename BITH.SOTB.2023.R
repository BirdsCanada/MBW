# For Trend estimation: https://besjournals.onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2009.01724.x
# and comment here (https://bookdown.org/adam_smith2/bbsbayes_intro_workshop/IndicesTrends.html)
# "...fitting a log-linear slope to the series of all annual indices between the two end-points 
# (e.g., all 11 years in a 10-year trend from 2008-2018). The slope of this line could be expressed 
# as an average annual percent change across the time-period of interest." 
# also some good ideas on that page about plotting trends on the map for each region


# 2010-2023 State of the Birds Report
# BITH Bicknell's Thrush
# Fitting same abundance-side model for all species
# file is sort by year, then region (to facilitate regional trends), and route.stop
# I removed the regional intercept for three reasons
# 1. It's based on arbitrary political boundaries (and we already have Latitude and Elevation)
# 2. DIC suggests that it's not a better model than a simpler one with just annual intercepts
# 3. Although removing the regional intercept changes the Abundance-Latitude relationship, the trends don't seem to change much
# I also considered Longitude, but Lat*Lon is highly (0.72) correlated

#Libraries
library(jagsUI)
library(beepr)
library(wesanderson)
library(colorspace)

rm(list=ls())
par(mfrow=c(1,1),mar=c(5,5,3,2)+0.1,cex.axis=1.5,cex.lab=1.5)
setwd('c:/Users/JasonHill/Documents/R/MBW_SOTB/2023/BITH')

#Import data
dff<-read.csv('BITH.SOTB2023.Kery.Format.csv',header=T,na.strings = "NA",stringsAsFactors = 1);str(dff)
df<-read.table('BITH.SOTB2023.Counts.txt',header=T,stringsAsFactors = 1);str(df) # BITH counts
dfcov<-read.table('AnySpecies.SOTB2023.CountTime.txt',header=T);str(dfcov) # sampling covariate count time (CTSR1 through CTSR4), expressed as minutes before or after local sunrise

# Put counts into blank array
ydata<-array(NA,dim=c(max(dff$Route.Stop2), 4, max(dff$Year2)));str(ydata) # 791 sampling locations, 4 reps (i.e., 5-minute point counts)/day, 14 years
for (k in 1:max(dff$Year2)) {
  sel.rows <- df$Year2==k
  ydata[,,k] <- as.numeric(as.matrix(df)[sel.rows,3:6]) #grab columns 3-6
}; head(ydata) # looks correct

#############################################################################################################################
# Now bring in covariate data
#############################################################################################################################

# Elevation doesn't change so we only need 1 value for each of the sampling stations
ELEV<-dff$Elevation[1:max(dff$Route.Stop2)] ;hist(ELEV)
ELEV.sc<-(ELEV-mean(ELEV))/sd(ELEV) ;hist(ELEV.sc) # center and scale elevation covariate to help JAGS run well
ELEV2<-ELEV.sc^2 ;hist(ELEV2) #quadratic term
ELEV3<-ELEV.sc^3 ;hist(ELEV3) #3rd term

# we'll create a matrix of elevation data for plotting/predicting
# since we have an intercept for region we'll end up arbitrarily plotting the elevation for some region
# we'll just pick Vermont, so we should constrain elevation values to Vermont when we plot it perhaps
ELEV.data<-matrix(NA,nrow=100,ncol=4)
ELEV.data[,1]<-seq(min(ELEV),max(ELEV),,100)
ELEV.data[,2]<-(ELEV.data[,1]-mean(ELEV))/sd(ELEV) # elevation data scaled and centered
ELEV.data[,3]<-ELEV.data[,2]^2 # quadratic elevation term
ELEV.data[,4]<-ELEV.data[,2]^3 #3rd term

# Latitude
Latitude.sc<-(dff$Latitude[1:max(dff$Route.Stop2)]-mean(dff$Latitude))/sd(dff$Latitude) ; hist(Latitude.sc)
Latitude2<-Latitude.sc^2 ;hist(Latitude2) #create quadratic term
Latitude3<-Latitude.sc^3 ;hist(Latitude3) #create quadratic term

Lat.data<-matrix(NA,nrow=100,ncol=4) #create matrix to make plotting predictions easier later
Lat.data[,1]<-seq(min(dff$Latitude),max(dff$Latitude),,100)
Lat.data[,2]<-(Lat.data[,1]-mean(dff$Latitude))/sd(dff$Latitude)
Lat.data[,3]<-Lat.data[,2]^2
Lat.data[,4]<-Lat.data[,2]^3
head(Lat.data);hist(Lat.data[,2])

#################################### Count Time for each survey
# Create a blank array to hold our sampling covariate count time
CountTime<-array(NA,dim=c(max(dfcov$Route.Stop2),4,max(dfcov$Year2)))    #(this is the same format as our count data)
str(CountTime)
# now fill in the CountTime array
for (k in 1:max(dfcov$Year2)) {
  sel.rows <- dfcov$Year2==k
  CountTime[,,k] <- as.numeric(as.matrix(dfcov)[sel.rows,3:6])
}
summary(CountTime);hist(CountTime)#lots of missing values, mean is 54.37
sd(CountTime,na.rm=1) #62.04143
head(CountTime[,,1:2])
# note the missing values (NA); that will cause all kinds of problems in JAGS
# 3 ways to do this:.
#1) replace all missing values with 0 or the mean
#2) replace all missing values by sampling from the other non-NA values already in the dataset; this would keep the shape of the distribution approximately the same
#3) replace all missing values for first point count period w/ mean, then add 5 minutes to each point count period
# I've messed around and none of these approaches changes anything noticeable in the final model results,
# and the model doesn't run any faster with all missing values set to zero (bummer)

#option 1
#CountTime[is.na(CountTime)]<-mean(CountTime,na.rm=1) # or just <-0
#CT.sc<-(CountTime-mean(CountTime,na.rm=1))/sd(CountTime,na.rm=1);hist(CT.sc)
#CT2<-CT.sc^2;hist(CT2)
#summary(CT.sc);summary(CT2)

#option 2
#CTvals<-CountTime[!is.na(CountTime)];summary(CTvals) #create vector of non-NA values to sample from
#CountTime[is.na(CountTime)]<-sample(CTvals,sum(is.na(CountTime)),replace=1)
#length(CountTime[is.na(CountTime)]) # shouldn't be anymore NAs left

#option 3, the one we're going with for now
mean(CountTime[,1,],na.rm=TRUE) #46.87108
mean(c(39,44,49,54)) # close enough
for (i in 1:4){
  CountTime[,i,][is.na(CountTime[,i,])]<-39+i*5 #add 5 minutes each period
}
CT.sc<-(CountTime-mean(CountTime))/sd(CountTime);hist(CT.sc)
CT2<-CT.sc^2;hist(CT2);summary(CT2) # No NAs
CT3<-CT.sc^3;hist(CT3);summary(CT3) # No NAs

# the CT.data is just for predictions, and we're only interested in showing predictions for data collected from 45 minutes before local sunrise
# to ~200 minutes after sunrise [only a tiny portion of the data comes from after 8:20 am]
CT.data<-matrix(NA,nrow=100,ncol=4)
CT.data[,1]<-seq(-45,200,,100)
CT.data[,2]<-(CT.data[,1]-mean(CountTime))/sd(CountTime) # must center/scale our prediction data by the same values as our real data
CT.data[,3]<-CT.data[,2]^2
CT.data[,4]<-CT.data[,2]^3
(CT.data)

##################### YEAR
# the count data are arranged in an array, by year, so you don't need something like "alpha.year[YEAR[i]]"
# you just run a loop for k in 1:14 years
# but for ease of plotting we'll create a year matrix for use with the regional trends
YEAR.data<-matrix(NA,nrow=100,ncol=2)
YEAR.data[,1]<-seq(2010,2023,,100)
YEAR.data[,2]<-YEAR.data[,1]-2016 # centered

##################### WIND force
str(dff$WindForceID);hist(dff$WindForceID);median(dff$WindForceID,na.rm = TRUE);head(dff$WindForceID);table(dff$WindForceID)
#lots of missing values (from non-survey events (of course), but observers did not always record those values))
# must fill in; filling in with median value (1, in this case) would be easy....like the next line of code
#WIND<-dff$WindForceID[is.na(dff$WindForceID)]<-1 #if dff$WindForceID is NA then replace with 1

# but instead, we'll sample from the actual non-NA Wind value distribution like so...
#option 2
WIND.Values<-dff$WindForceID[!is.na(dff$WindForceID)];summary(WIND.Values);hist(WIND.Values);table(WIND.Values) #create vector of non-NA values to sample from
windy<-dff$WindForceID
windy[is.na(windy)]<-sample(WIND.Values,sum(is.na(windy)),replace=1)
length(windy[is.na(windy)]);hist(windy);table(windy) # shouldn't be anymore NAs left

#looks good--now we need to make it into a multi-year array (like CountTime) since wind values change at each sampling location each year

# Let's create a blank array to hold our sampling covariate WIND
WIND<-array(NA,dim=c(max(dfcov$Route.Stop2),max(dfcov$Year2)));str(WIND)    # similar to CountTime but only two dimensions--wind value remains the same for all 4 point counts at a sampling station in a year

# now fill in the WIND array
for (k in 1:max(dfcov$Year2)) {
  WIND[,k] <- windy[(1+((k-1)*max(dfcov$Route.Stop2))):(max(dfcov$Route.Stop2)*k)]
}

#make into a standardized covariate
WIND.sc<-(WIND-mean(WIND))/sd(WIND);hist(WIND.sc)
WIND2<-WIND.sc^2;hist(WIND2);summary(WIND2)
WIND3<-WIND.sc^3;hist(WIND3);summary(WIND3)

# let's make an array for predictions
WIND.data<-matrix(NA,nrow=100,ncol=4)
WIND.data[,1]<-seq(0,5,,100) #100 values between 0 and 5
WIND.data[,2]<-(WIND.data[,1]-mean(WIND))/sd(WIND) # must center/scale our prediction data by the same values as our real data
WIND.data[,3]<-WIND.data[,2]^2
WIND.data[,4]<-WIND.data[,2]^3
(WIND.data)

##################### Noise?
# The Noise value is more tricky
# it mostly consists of missing values--b/c we just started collecting it in 2019
# we'll ignore for now, but most likely Noise values aren't random--some sites are more likely to be windy every year
# I think we'd use actual wind values (when they were recorded), and then for missing values
# we'd estimate the missing value based on other values from that same sampling point in other years.
# But what if there hasn't ever been a noise value for that sampling station (e.g., some routes were retired before we started collecting those noise ratings)
# In that case, we could estimate Noise values based on a model including historical weather data and/or elevation/lat/long?
# I think we leave that for another day

cor.test(dff$WindForceID,dff$Noise,use='complete.obs')
#correlation of .54 between noise and wind force....hmmmm
# wind speed isn't related to noise from water though

##################### June Day
str(dff$JuneDay);hist(dff$JuneDay);median(dff$JuneDay,na.rm = TRUE);head(dff$JuneDay)
#lots of missing values (from non-survey events
# must fill in; filling in with median value (17) is easy

JD<-dff$JuneDay
JD[is.na(JD)]<-17;summary(JD)

#looks good--now we need to make it into a multi-year array (like WIND) since the survey date values change at each sampling location each year

# Let's create a blank array to hold our sampling covariate JuneDay
JuneDay<-array(NA,dim=c(max(dfcov$Route.Stop2),max(dfcov$Year2)));str(JuneDay)    # similar to CountTime but only two dimensions--JuneDay value remains the same for all 4 point counts at a sampling station in a year

# now fill in the JuneDay array
for (k in 1:max(dfcov$Year2)) {
  JuneDay[,k] <- JD[(1+((k-1)*max(dfcov$Route.Stop2))):(max(dfcov$Route.Stop2)*k)]
}

#make into a standardized covariate
JD.sc<-(JuneDay-mean(JuneDay))/sd(JuneDay);hist(JD.sc)
JD2<-JD.sc^2;hist(JD2);summary(JD2)
JD3<-JD.sc^3;hist(JD3);summary(JD3)

# let's make an array for predictions
JD.data<-matrix(NA,nrow=100,ncol=4)
JD.data[,1]<-seq(1,30,,100) #100 values between 1 and 30
JD.data[,2]<-(JD.data[,1]-mean(JuneDay))/sd(JuneDay) # must center/scale our prediction data by the same values as our real data
JD.data[,3]<-JD.data[,2]^2
JD.data[,4]<-JD.data[,2]^3
(JD.data)

#############################################################################################################################
# Now we can specify the model. The code below doesn't show this, but I start with a crazy basic model and then add in 1 additional parameter at a time
# Kery and Royle refer to this as "sneaking up' on the model
#############################################################################################################################
sink("BITH_SOTB2023.txt")
cat("
    model {
    
# Priors for detection and abundance intercepts 
# intercepts for years 1 through 14, to allow mean detection rates to vary by year
# you would expect them to differ by year just based on weather, observer differences, etc.
for (k in 1:Y) {
    alpha.det.year[k] ~ dnorm(0,0.01)
    alpha.lam.year[k] ~ dnorm(0,tau.year)                   # random year effect/intercept
    }
    tau.year<-1/(sd.year*sd.year)
    sd.year~dunif(0,1)                                      # sd for random year effect on abundance

  for (r in 1:5) { # 5 regions (NY_Catskills, NY_Adriondacks, VT, NH, ME based on management units)
  alpha.lam.trend.region[r] ~ dnorm(mu.int, tau.int)	   # Random intercepts
  beta.lam.trend.region[r] ~ dnorm(mu.slope, tau.slope)  # Random slopes
}

mu.int ~ dnorm(0, 0.01)		# Mean hyperparameter for random intercepts
tau.int <- 1 / (sigma.int * sigma.int)
sigma.int ~ dunif(0, 10)		# SD hyperparameter for random intercepts

mu.slope ~ dnorm(0, 0.01)		# Mean hyperparameter for random slopes
tau.slope <- 1 / (sigma.slope * sigma.slope)
sigma.slope ~ dunif(0, 5)		# SD hyperparameter for slopes 

# Variance terms for site random effect
    tau.lam <- 1 / (sd.lam * sd.lam)
    sd.lam ~ dunif(0, 3)     # SD of random site effect on abundance 

# Prior for detection parameters
    beta.det.CT ~ dlogis(0,1)                               # Count.time
    beta.det.CT2 ~ dlogis(0,1)                              # Count.time2
    beta.det.JD ~ dlogis(0,1)                               # JuneDay
    beta.det.JD2 ~ dlogis(0,1)                              # JuneDay2
    beta.det.WIND ~ dlogis(0,1)                             # Wind
    beta.det.WIND2 ~ dlogis(0,1)                            # Wind2
    
# Priors for abundance (lambda) parameters
    beta.lam.lat ~ dlogis(0,1)                              # weakly informative prior fo latitude effect
    beta.lam.lat2 ~ dlogis(0,1)                              # weakly informative prior fo latitude2 effect
    beta.lam.elev ~ dlogis(0,1)                             # weakly-informative prior for elevation effect
    beta.lam.elev2 ~ dlogis(0,1)                             # weakly-informative prior for elevation2 effect

# Likelihood, ecological model for true abundance
    for (i in 1:R) {                                        # loop over R (791) sites
    eps[i] ~ dnorm(0, tau.lam)                             # random site effect, abundance noise
    for (k in 1:Y) {                                        # loop over k (14) years
    N[i,k] ~ dpois(lambda[i,k])                             # Abundance [site i, year k]

# assemble the abundance component of the model
    log(lambda[i,k]) <- alpha.lam.trend.region[REGION[i]] + beta.lam.trend.region[REGION[i]]*year[k] +
    beta.lam.elev*COV.elev[i] + beta.lam.elev2*COV.elev2[i] + 
    beta.lam.lat*COV.lat[i] + beta.lam.lat2*COV.lat2[i] + 
    alpha.lam.year[k] + eps[i]
    
# Observational model for replicated counts
    for (j in 1:reps) {                                     # loop over temporal replicates (4 point counts per sampling station per year)
    ydata[i,j,k] ~ dbin(p[i,j,k],N[i,k])                    # Detection & abundance: indexed by site, survey, and year
    p[i,j,k] <- exp(lp[i,j,k])/(1+exp(lp[i,j,k]))           # define logit function here
# assemble the detecion components
    lp[i,j,k] <-alpha.det.year[k] + beta.det.CT*COV.CT[i,j,k] + beta.det.CT2*COV.CT2[i,j,k] + 
    beta.det.JD*COV.JD[i,k] + beta.det.JD2*COV.JD2[i,k] +
    beta.det.WIND*COV.WIND[i,k] + beta.det.WIND2*COV.WIND2[i,k]
    
# Assess model fit using Chi-square discrepancy
  # Compute fit statistic E for observed data
#    eval[i,j,k] <- p[i,j,k] * N[i,k] # Expected values
#    E[i,j,k]    <- pow((ydata[i,j,k] - eval[i,j,k]),2)/(eval[i,j,k] + 0.005) # add 0.005 to prevent any chance of division by zero
    
# Generate replicate data and compute fit stats for them
#    y.new[i,j,k] ~ dbin(p[i,j,k],N[i,k])
#    E.new[i,j,k] <- pow((y.new[i,j,k] - eval[i,j,k]),2) / (eval[i,j,k] + 0.005) # add 0.005 to prevent division by zero
    }}} #j,k,i
    
# Derived and other quantities
    for (k in 1:Y) {
#    # N gets estimated for every site in every year, including for sites that didn't have point counts performed at them
    totalN[k] <- sum(N[,k])                                 # Total pop size across all sites within each year...might call this the pop size surrounding sampling sites
    Maine.N[k]<- sum(N[1:187,k])
    NewHampshire.N[k]<- sum(N[188:466,k])
    NewYorkA.N[k]<- sum(N[467:585,k])
    NewYorkC.N[k]<- sum(N[586:626,k])
    NewYork.N[k]<- sum(N[467:626,k])
    Vermont.N[k]<- sum(N[627:791,k])
    }

#    fit <- sum(E[,,])
#    fit.new <- sum(E.new[,,])

# use the mean detection and abundance betas in these plots, so the plot shows elevation and count time relationships in a typical year
    for (b in 1:100) {
    lam.pred.elev[b]<-exp(beta.lam.elev * elev.pred[b] + beta.lam.elev2 * elev2.pred[b])
    lam.pred.lat[b]<-exp(beta.lam.lat * lat.pred[b] + beta.lam.lat2 * lat2.pred[b])
    lam.pred.NE.trend[b]<-exp(mu.int + mu.slope * year.pred[b])
    lam.pred.ME.trend[b]<-exp(alpha.lam.trend.region[1] + beta.lam.trend.region[1] * year.pred[b] + beta.lam.elev * -0.7466614 + beta.lam.elev2 * 0.5575 + beta.lam.lat * 1.1636 + beta.lam.lat2 * 1.354098)
    lam.pred.NH.trend[b]<-exp(alpha.lam.trend.region[2] + beta.lam.trend.region[2] * year.pred[b] + beta.lam.elev * 0.0789 + beta.lam.elev2 * 0.00623 + beta.lam.lat * -0.1 + beta.lam.lat2 * 0.01)
    lam.pred.NYA.trend[b]<-exp(alpha.lam.trend.region[3] + beta.lam.trend.region[3] * year.pred[b] + beta.lam.elev * 0.605173 + beta.lam.elev2 * .3662 + beta.lam.lat * -0.14658 + beta.lam.lat2 * 0.021487)
    lam.pred.NYC.trend[b]<-exp(alpha.lam.trend.region[4] + beta.lam.trend.region[4] * year.pred[b] + beta.lam.elev * 0.8931 + beta.lam.elev2 * .7976 + beta.lam.lat * -2.8762 + beta.lam.lat2 * 8.2724)
    lam.pred.VT.trend[b]<-exp(alpha.lam.trend.region[5] + beta.lam.trend.region[5] * year.pred[b] + beta.lam.elev * 0.0544 + beta.lam.elev2 * 0.0029 + beta.lam.lat * -0.3293 + beta.lam.lat2 * 0.1084)
    logit(p.pred.CT[b])<-alpha.det.year[9] + beta.det.CT * CT.pred[b] + beta.det.CT2 * CT2.pred[b]
    logit(p.pred.JD[b])<-alpha.det.year[9] + beta.det.JD * JD.pred[b] + beta.det.JD2 * JD2.pred[b]
    logit(p.pred.WIND[b])<-alpha.det.year[9] + beta.det.WIND * WIND.pred[b] + beta.det.WIND2 * WIND2.pred[b]
    }
  }
    ",fill=TRUE)
sink()

############################################################################################
# Bundle the data
R = nrow(ydata);R
reps = ncol(ydata);reps
Y = max(df$Year2);Y

COV.elev=ELEV.sc ;hist(COV.elev)
COV.elev2=ELEV2 ;hist(COV.elev2)
COV.CT<-CT.sc; hist(COV.CT)
COV.CT2<-CT2; hist(COV.CT2)
COV.lat<-Latitude.sc;hist(COV.lat)
COV.lat2<-Latitude2;hist(COV.lat2)
COV.JD<-JD.sc;hist(COV.JD)
COV.JD2<-JD2;hist(COV.JD2)
COV.WIND<-WIND.sc;hist(COV.WIND)
COV.WIND2<-WIND2;hist(COV.WIND2)
REGION<-as.factor(dff$Region[1:791]);summary(REGION)

# Create prediction data sets to facilitate plotting
elev.pred<-ELEV.data[,2] ;hist(elev.pred)
elev2.pred<-ELEV.data[,3] ;hist(elev2.pred)
lat.pred<-Lat.data[,2] ;hist(lat.pred)
lat2.pred<-Lat.data[,3] ;hist(lat2.pred)
CT.pred<-CT.data[,2] ;hist(CT.pred)
CT2.pred<-CT.data[,3] ;hist(CT2.pred)
WIND.pred<-WIND.data[,2] ;hist(WIND.pred)
WIND2.pred<-WIND.data[,3] ;hist(WIND2.pred)
JD.pred<-JD.data[,2] ;hist(JD.pred)
JD2.pred<-JD.data[,3] ;hist(JD2.pred)
year<-(1:14)-7
year.pred<-YEAR.data[,2];hist(year.pred)

#Assemble the data for JAGS
JAGS.data <- list(ydata=ydata,R=R,reps=reps,Y=Y,year.pred=year.pred,
                  COV.elev=COV.elev,COV.elev2=COV.elev2,
                  COV.CT=COV.CT,COV.CT2=COV.CT2,
                  COV.JD=COV.JD,COV.JD2=COV.JD2,
                  COV.WIND=COV.WIND,COV.WIND2=COV.WIND2,
                  year=year,REGION=as.numeric(REGION),
                  COV.lat=COV.lat,COV.lat2=COV.lat2,
                  lat.pred=lat.pred,lat2.pred=lat2.pred,
                  CT.pred=CT.pred,CT2.pred=CT2.pred,
                  JD.pred=JD.pred,JD2.pred=JD2.pred,
                  WIND.pred=WIND.pred,WIND2.pred=WIND2.pred,
                  elev.pred=elev.pred,elev2.pred=elev2.pred)

# Generate initial values
Nst <- apply(ydata,c(1,3),max) + 1;  str(Nst)     # add 1 to max count as estimate of site-specific N
Nst[is.na(Nst)] <- 1                              # supply initial values of 1 for each missing count, or JAGS will throw up violently

inits <- function(){list(N=Nst,alpha.lam.year=rep(0,14),sd.year=0.1,
                         beta.det.CT=-0.1,beta.det.CT2=-0.1,
                         beta.det.WIND=-0.1,beta.det.WIND2=0,
                         beta.lam.elev=0.3,beta.lam.elev2=-.2,
                         beta.det.JD=-0.1,beta.det.JD2=0,
                         beta.lam.lat=0.16,beta.lam.lat2=0,
                         alpha.lam.trend.region=rep(-.5,5),
                         beta.lam.trend.region=rep(-.02,5),
                         mu.int=-.5,sigma.int=.5,
                         mu.slope=-0.02,sigma.slope=.05,sd.lam=1
)}

# List the parameters that we want to monitor
params <- c("totalN","sd.year","sd.lam","sigma.slope","sigma.int",
            "alpha.lam.trend.region","mu.int","beta.lam.trend.region","mu.slope",
            "Maine.N","NewHampshire.N","NewYorkA.N","NewYorkC.N","NewYork.N","Vermont.N",
            "lam.pred.ME.trend","lam.pred.NH.trend",
            "lam.pred.NYA.trend","lam.pred.NYC.trend","lam.pred.VT.trend","lam.pred.NE.trend",
            "beta.lam.elev","beta.lam.elev2",
            "beta.lam.lat","beta.lam.lat2",
            "beta.det.CT","beta.det.CT2",
            "beta.det.WIND","beta.det.WIND2",
            "beta.det.JD","beta.det.JD2",
            "alpha.det.year","alpha.lam.year",
#           "fit","fit.new",
            "p.pred.JD","p.pred.CT","p.pred.WIND","lam.pred.elev","lam.pred.lat","lam.pred.WIND")


# Specify the MCMC Settings
ni<-50
nt<-1
nb<-00
nc<-3

ni<-10500
nt<-10
nb<-500
nc<-3

((ni-nb)/nt)    # of iterations saved per chain
((ni-nb)/nt)*nc # total number of iterations saved

# This model is fit, and ready to run
ptm_start <- Sys.time() #an easy way to record computational time
outSOTBBITH<-jags(data = JAGS.data,inits=inits,parameters.to.save=params,parallel=TRUE,
                  model.file = "BITH_SOTB2023.txt",n.chains=nc, n.thin=nt, n.iter=ni, n.burnin=nb);ptm_end<-Sys.time() - ptm_start; ptm_end # 
ptm_start <- Sys.time();beep(sound=8)

pp.check(outSOTBBITH,observed="fit",simulated = "fit.new",pch=16,col="#fdae6b")
write.csv(outSOTBBITH$summary,"Summary_TrendModel.SOTBBITH.csv")
write.csv(outSOTBBITH$sims.list,"Sims_TrendModel.SOTBBITH.csv")
save.image(file="TrendModel.BITH.JAGS.RData") 
#load(file="TrendModel.BITH.JAGS.RData") # in case you ever need to reload the .RData
####################################################################################################################################################
# Results plots, draw samples
par(mfrow=c(1,1),mar=c(5,7,3,2)+0.1,cex.axis=1.1,cex.lab=1)
mcmc.sample<-((ni-nb)/nt)*nc
sub.set<-sort(sample(1:mcmc.sample,size=1500))
####################################################################################################################################################
# Abundance plot
# Should we consider smoothing this line?
png("BITH.StudyAreaAbundance.png", width = 6*.8, height = 4*.8,units="in",res=500)
par(mfrow=c(1,1),mar=c(3,7,1,1)+0.1,cex.axis=1.1,cex.lab=1)
plot(1:14,outSOTBBITH$mean$totalN,xlab="",ylab="",las=1,ylim=c(200,600),type = "l",frame.plot = FALSE,pch=16,axes = 0,
     main="",col= wes_palette("Royal1")[4],lwd=3,cex.axis=.8)
axis(1, at=1:14, labels=c(2010:2023),cex.axis=.8)
axis(2, las=1,cex.axis=.8)
mtext(side=2,text="Bicknell's Thrush",line=3.75,cex=.9,col=darken(wes_palette("Royal1")[4],.3))
mtext(side=2,text="local population size", line=2.5, cex=.9,col=darken(wes_palette("Royal1")[4],.3))
segments(1:14,outSOTBBITH$q2.5$totalN,1:14,outSOTBBITH$q97.5$totalN,col=lighten(wes_palette("Royal1")[4],.3),lwd=3)
lines(1:14,outSOTBBITH$mean$totalN,col=wes_palette("Royal1")[4],lwd=3)
dev.off()


png("BITH.StudyAreaAbundance.png", width = 6*.8, height = 4*.8,units="in",res=700)
par(mfrow=c(1,1),mar=c(2,6,1,1)+0.1,cex.axis=1.1,cex.lab=1)
plot(1:14,outSOTBBITH$mean$totalN,xlab="",ylab="",las=1,ylim=c(200,600),type = "n",frame.plot = FALSE,pch=16,axes = 0,
     main="",col= wes_palette("Royal1")[4],lwd=3,cex.axis=.5)
axis(1, at=1:14, labels=FALSE,cex.axis=.5,cex=.5,tck=-0.02)
mtext(side=1,text=c("2010","2014","(Survey year)","2019","2023"),at=c(1,5,7.5,10,14),line=0.2,cex=.5)
axis(2, las=1,cex.axis=.8)
mtext(side=2,text="Bicknell's Thrush",line=3.75,cex=.9,col=darken(wes_palette("Rushmore1")[3],.3))
mtext(side=2,text="local population size", line=2.5, cex=.9,col=darken(wes_palette("Rushmore1")[3],.3))
#for (i in sub.set){
#  lines(1:14,outSOTBBITH$sims.list$totalN[i,],
#         lwd=1,col=adjustcolor(col=wes_palette("Rushmore1")[3],alpha.f=.05))}
for (i in sub.set){
  lines(smooth.spline(1:14,outSOTBBITH$sims.list$totalN[i,],spar=.2),
         lwd=1,col=adjustcolor(col=wes_palette("Rushmore1")[3],alpha.f=.05))}

#lines(1:14,outSOTBBITH$mean$totalN,col=wes_palette("Rushmore1")[4],lwd=3)
# or
lines(smooth.spline(1:14, outSOTBBITH$mean$totalN, spar=0.2),lwd=3,col=wes_palette("Rushmore1")[4])
dev.off()
#






plot(outSOTBBITH$mean$lam.pred.ME.trend,pch=16,col="PURPLE",ylim=c(0,6),axes=0,ylab="",type="l",lwd=5)
lines(outSOTBBITH$mean$lam.pred.NH.trend,pch=16,col="BLUE",lwd=5)
lines(outSOTBBITH$mean$lam.pred.NYC.trend,pch=16,col="GREEN",lwd=5)
lines(outSOTBBITH$mean$lam.pred.NYA.trend,pch=16,col="ORANGE",lwd=5)
lines(outSOTBBITH$mean$lam.pred.VT.trend,pch=16,col="RED",lwd=5)
lines(outSOTBBITH$mean$lam.pred.NE.trend,pch=16,col="GRAY",lwd=5)
axis(1,at=c(1,17,32,47,62,78,93,100),labels=c(2010,2012,2014,2016,2018,2020,2022,""))
axis(2,las=1)
mtext(side=2,text="Bicknell's Thrush predicated abundance",line=3.75, cex=.85)
mtext(side=2,text="at a typical sampling station",line=2.75, cex=.85)



####################################################################################################################################################
# Count Time
#jpeg("BITH.Count.time.jpeg", width = 6, height = 4,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,3,2)+0.1,cex.axis=1.1,cex.lab=1)
plot(CT.pred,outSOTBBITH$mean$p.pred.CT,xlab="",ylab="",ylim=c(0.1,1),frame=F,type="l",las=1,axes=0,cex.lab=1.05,xlim=c(-2,max(CT.pred)),col=wes_palette("Rushmore1")[3],lwd=3)
sub.set<-sort(sample(1:mcmc.sample,size=200))
for (i in sub.set){
  points(CT.data[,2],jitter(outSOTBBITH$sims.list$p.pred.CT[i,]),
         type="l",lwd=1,col=adjustcolor(col=wes_palette("Rushmore1")[3],alpha.f=.05))}
sub.set<-sort(sample(1:mcmc.sample,size=200))
for (i in sub.set){
  points(CT.data[,2],jitter(outSOTBBITH$sims.list$p.pred.CT[i,]),
         type="l",lwd=1,col=adjustcolor(col=wes_palette("Rushmore1")[3],alpha.f=.05))}
lines(CT.pred,outSOTBBITH$mean$p.pred.CT,lty=1,lwd=4,col=wes_palette("Rushmore1")[3])
axis(2,las=1,cex.axis=.85)
axis(1,labels=c(-45,15,75,145,195),at=c(-1.91,-0.748,0.42,1.59,2.75),cex.axis=.85)
mtext(side=2,text="Bicknell's Thrush",line=4.5,cex=.85,col=wes_palette("Rushmore1")[3])
mtext(side=2,text=expression(paste("detection probability  ( ", italic(p),")")),line=3,cex=.85,,col=wes_palette("Rushmore1")[3])
mtext(side=1,text="Point count start time (minutes relative to local sunrise)",line=2.5,cex=.85,,col=wes_palette("Rushmore1")[3])
#dev.off()

# or for a more traditional plot
jpeg("BITH.Count.time.traditional.jpeg", width = 6, height = 4,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,3,2)+0.1,cex.axis=1.1,cex.lab=1)
plot(CT.pred,outSOTBBITH$mean$p.pred.CT,xlab="",ylab="", ylim=c(0,1),frame=F,type="l",las=1,axes=F)
polygon(c(CT.pred,rev(CT.pred)),c(outSOTBBITH$q2.5$p.pred.CT,rev(outSOTBBITH$q97.5$p.pred.CT)),col=lighten(wes_palette("Rushmore1")[3],.4),border=F)
lines(CT.pred,outSOTBBITH$mean$p.pred.CT,lty=1,lwd=3,col=wes_palette("Rushmore1")[3])
axis(2,las=1,cex.axis=0.85)
axis(1,labels=c(-45,15,75,135,195),at=c(-1.91,-0.748,0.42,1.59,2.75),cex.axis=.85)
mtext(side=2,text="Bicknell's Thrush",line=4.5,cex=.85,col=wes_palette("Rushmore1")[3])
mtext(side=2,text=expression(paste("detection probability ( ", italic(p),")")),line=3,cex=.85,col=wes_palette("Rushmore1")[3])
mtext(side=1,text="Point count start time (minutes relative to local sunrise)",line=2.5,cex=.85)
dev.off()

#################################################################################################################################
# Elevation relationship
png("BITH.Elevation.png", width = 6*.8, height = 4*.8,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,1,1)+0.1,cex.axis=1.1,cex.lab=1)
plot(elev.pred,outSOTBBITH$mean$lam.pred.elev,xlab="",ylab="",ylim=c(0,1.5),frame=F,type="l",las=1,axes=0,col=wes_palette("FantasticFox1")[3],lwd=3)
for (i in sub.set){
  points(elev.pred,jitter(outSOTBBITH$sims.list$lam.pred.elev[i,]),
         type="l",lwd=1,col=adjustcolor(col=wes_palette("FantasticFox1")[3],alpha.f=.1))}
sub.set<-sort(sample(1:mcmc.sample,size=200))
for (i in sub.set){
  points(elev.pred,jitter(outSOTBBITH$sims.list$lam.pred.elev[i,]),
         type="l",lwd=1,col=adjustcolor(col=wes_palette("FantasticFox1")[3],alpha.f=.1))}
lines(elev.pred,outSOTBBITH$mean$lam.pred.elev,lty=1,lwd=3,col=darken(wes_palette("FantasticFox1")[3],.5))
axis(2,las=1,cex.axis=.9,labels=c("0.0","0.5","1.0","1.5"),at=c(0,.5,1,1.5))
axis(1,labels=c(550,867,1184,1500),at=c(-3.02,-0.92,1.16,3.254),cex.axis=.9)
mtext(side=2,text="Bicknell's Thrush",line=4.1, cex=.75,col=darken(wes_palette("FantasticFox1")[3],.5))
mtext(side=2,text="predicted abundance at a",line=3.35,cex=.75,col=darken(wes_palette("FantasticFox1")[3],.5))
mtext(side=2,text="typical sampling station",line=2.5, cex=.75,col=darken(wes_palette("FantasticFox1")[3],.5))
mtext(side=1,text="Elevation (m)",line=2,cex=.9)
dev.off()

# or for a more traditional plot
jpeg("BITH.Elevation.traditional.jpeg", width = 6, height = 4,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,1,1)+0.1,cex.axis=1.1,cex.lab=1)
plot(elev.pred,outSOTBBITH$mean$lam.pred.elev,xlab="",ylab="", ylim=c(0,12),frame=F,type="l",las=1,axes=F)
polygon(c(elev.pred,rev(elev.pred)),c(outSOTBBITH$q2.5$lam.pred.elev,rev(outSOTBBITH$q97.5$lam.pred.elev)),col=lighten(wes_palette("FantasticFox1")[3],.4),border=F)
lines(elev.pred,outSOTBBITH$mean$lam.pred.elev,lty=1,lwd=3,col=wes_palette("FantasticFox1")[3])
axis(2,las=1,cex.axis=.85)
axis(1,labels=c(550,740,930,1120,1310,1500),at=c(-3.02,-1.766,-0.51,0.74,1.99,3.254),cex.axis=.85)
mtext(side=2,text="Bicknell's Thrush",line=4.5,cex=.85,col=darken(wes_palette("FantasticFox1")[3],.4))
mtext(side=2,text="mean expected abundance",line=3,cex=.85,col=darken(wes_palette("FantasticFox1")[3],.4))
mtext(side=1,text="Elevation (m)",line=2.5,cex=.85)
dev.off()

#################################################################################################################################
# let's plot the Latitude relationship
png("BITH.Latitude.png", width = 6*.8, height = 4*.8,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,1,1)+0.1,cex.axis=1.1,cex.lab=1)
plot(lat.pred,outSOTBBITH$mean$lam.pred.lat,xlab="",ylab="",ylim=c(0,2),frame=F,type="l",las=1,axes=0,col=wes_palette("GrandBudapest2")[1],lwd=3)
for (i in sub.set){
  points(lat.pred,jitter(outSOTBBITH$sims.list$lam.pred.lat[i,]),
         type="l",lwd=1,col=adjustcolor(col=wes_palette("GrandBudapest2")[1],alpha.f=.1))}
sub.set<-sort(sample(1:mcmc.sample,size=200))
for (i in sub.set){
  points(lat.pred,jitter(outSOTBBITH$sims.list$lam.pred.lat[i,]),
         type="l",lwd=1,col=adjustcolor(col=wes_palette("GrandBudapest2")[1],alpha.f=.1))}
lines(lat.pred,outSOTBBITH$mean$lam.pred.lat,lty=1,lwd=3,col=darken(wes_palette("GrandBudapest2")[1],.4))
axis(2,las=1,at=c(0,1,2),labels=c("0.0","1.0","2.0"),cex.axis=.9)
axis(1,labels=c("42.0","43.0","44.0","45.0","46.0"),at=c(-3.15,-1.81,-0.48,0.85,2.19),cex.axis=.9)
mtext(side=2,text="Bicknell's Thrush",line=4.1,cex=.85,col=darken(wes_palette("GrandBudapest2")[1],.5))
mtext(side=2,text="predicted abundance at a",line=3.35,cex=.85,col=darken(wes_palette("GrandBudapest2")[1],.5))
mtext(side=2,text="typical sampling station",line=2.5,cex=.85,col=darken(wes_palette("GrandBudapest2")[1],.5))
mtext(side=1,text=expression(~Latitude~degree~N),line=2.0,cex=.85)
dev.off()

# or for a more traditional plot
jpeg("BITH.Latitude.traditional.jpeg", width = 6, height = 4,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,1,1)+0.1,cex.axis=1.1,cex.lab=1)
plot(lat.pred,outSOTBBITH$mean$lam.pred.lat,xlab="",ylab="", ylim=c(0,4),frame=F,type="l",las=1,axes=F,xlim=c(lat.pred[1],2.2))
polygon(c(lat.pred,rev(lat.pred)),c(outSOTBBITH$q2.5$lam.pred.lat,rev(outSOTBBITH$q97.5$lam.pred.lat)),col=lighten(wes_palette("GrandBudapest2")[1],.6),border=F)
lines(lat.pred,outSOTBBITH$mean$lam.pred.lat,lty=1,lwd=3,col=wes_palette("GrandBudapest2")[1])
axis(2,las=1,cex.axis=.85)
axis(1,labels=c(42.0,43.0,44.0,45.0,46.0),at=c(-3.15,-1.81,-0.48,0.85,2.19),cex.axis=.85)
mtext(side=2,text="Bicknell's Thrush",line=4.5,cex=.85,col=darken(wes_palette("GrandBudapest2")[1],.3))
mtext(side=2,text="mean expected abundance",line=3,cex=.85,col=darken(wes_palette("GrandBudapest2")[1],.3))
mtext(side=1,text=expression(~Latitude~degree~N),line=2.5,cex=.85,col="BLACK")
dev.off()

####################################################################################################################################################
# Wind Speed
#jpeg("BITH.WIND.jpeg", width = 6, height = 4,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,3,2)+0.1,cex.axis=1.1,cex.lab=1)
plot(WIND.pred,outSOTBBITH$mean$p.pred.WIND,xlab="",ylab="",ylim=c(0.1,1),frame=F,type="l",las=1,axes=0,cex.lab=1.05,xlim=c(-2,5.25),col="#31a354",lwd=3)
for (i in sub.set){
  points(WIND.data[,2],jitter(outSOTBBITH$sims.list$p.pred.WIND[i,]),
         type="l",lwd=1,col=adjustcolor(col="#78c679",alpha.f=.1))}
sub.set<-sort(sample(1:mcmc.sample,size=200))
for (i in sub.set){
  points(WIND.data[,2],jitter(outSOTBBITH$sims.list$p.pred.WIND[i,]),
         type="l",lwd=1,col=adjustcolor(col="#78c679",alpha.f=.1))}
lines(WIND.pred,outSOTBBITH$mean$p.pred.WIND,lty=1,lwd=4,col="#006837")
axis(2,las=1,cex.axis=.85)
axis(1,labels=c(-45,15,75,135,195),at=c(-1.91,-0.748,0.42,1.59,2.75),cex.axis=.85)mtext(side=2,text="Bicknell's Thrush",line=4.5,cex=.85)
mtext(side=2,text=expression(paste("detection probability  ( ", italic(p),")")),line=3,cex=.85)
mtext(side=1,text="Beaufort Wind Scale",line=2.5,cex=.85)
#dev.off()

# or for a more traditional plot
#jpeg("BITH.WindSpeed.traditional.jpeg", width = 6, height = 4,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,3,2)+0.1,cex.axis=1.1,cex.lab=1)
plot(WIND.pred,outSOTBBITH$mean$p.pred.WIND,xlab="",ylab="", ylim=c(0,1),frame=F,type="l",las=1,axes=F)
polygon(c(WIND.pred,rev(WIND.pred)),c(outSOTBBITH$q2.5$p.pred.WIND,rev(outSOTBBITH$q97.5$p.pred.WIND)),col="#78c679",border=F)
lines(WIND.pred,outSOTBBITH$mean$p.pred.WIND,lty=1,lwd=3,col="#006837")
axis(2,las=1,cex.axis=0.85)
axis(1,labels=c(0,1,2,3,4,5),at=c(-1.12,-0.376,0.367,1.111,1.855,2.61),cex.axis=.85)
mtext(side=2,text="Bicknell's Thrush",line=4.5,cex=.85)
mtext(side=2,text=expression(paste("detection probability ( ", italic(p),")")),line=3,cex=.85)
mtext(side=1,text="Beaufort Wind Scale",line=2.5,cex=.85)
dev.off()

#################################################################################################################################
#################################################################################################################################
# June Day
#jpeg("BITH.JuneDay.jpeg", width = 6, height = 4,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,3,2)+0.1,cex.axis=1.1,cex.lab=1)
plot(JD.pred,outSOTBBITH$mean$p.pred.JD,xlab="",ylab="",ylim=c(0.1,1),frame=F,type="l",las=1,axes=0,cex.lab=1.05,xlim=c(-2,5.25),col="#31a354",lwd=3)
for (i in sub.set){
  points(JD.data[,2],jitter(outSOTBBITH$sims.list$p.pred.JD[i,]),
         type="l",lwd=1,col=adjustcolor(col="#78c679",alpha.f=.1))}
sub.set<-sort(sample(1:mcmc.sample,size=200))
for (i in sub.set){
  points(JD.data[,2],jitter(outSOTBBITH$sims.list$p.pred.JD[i,]),
         type="l",lwd=1,col=adjustcolor(col="#78c679",alpha.f=.1))}
lines(JD.pred,outSOTBBITH$mean$p.pred.JD,lty=1,lwd=4,col="#006837")
axis(2,las=1,cex.axis=.85)
axis(1,labels=c(1,10,20,30),at=c(-2.135,-0.966,0.332,1.63),cex.axis=.85)
mtext(side=2,text=expression(paste("detection probability  ( ", italic(p),")")),line=3,cex=.85)
mtext(side=1,text="June Day (1 = June 1st)",line=2.5,cex=.85)
#dev.off()

# or for a more traditional plot
jpeg("BITH.JuneDay.traditional.jpeg", width = 6, height = 4,units="in",res=500)
par(mfrow=c(1,1),mar=c(5,7,3,2)+0.1,cex.axis=1.1,cex.lab=1)
plot(JD.pred,outSOTBBITH$mean$p.pred.JD,xlab="",ylab="", ylim=c(0,1),frame=F,type="l",las=1,axes=F)
polygon(c(JD.pred,rev(JD.pred)),c(outSOTBBITH$q2.5$p.pred.JD,rev(outSOTBBITH$q97.5$p.pred.JD)),col="#78c679",border=F)
lines(JD.pred,outSOTBBITH$mean$p.pred.JD,lty=1,lwd=3,col="#006837")
axis(2,las=1,cex.axis=0.85)
axis(1,labels=c(1,10,20,30),at=c(-2.135,-0.966,0.332,1.63),cex.axis=.85)
mtext(side=2,text="Bicknell's Thrush",line=4.5,cex=.85)
mtext(side=2,text=expression(paste("detection probability ( ", italic(p),")")),line=3,cex=.85)
mtext(side=1,text="June Day (1 = June 1st)",line=2.5,cex=.85)
dev.off()

#################################################################################################################################