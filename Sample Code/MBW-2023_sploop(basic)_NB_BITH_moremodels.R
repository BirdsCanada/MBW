
### species specific loop function
### This is for ALL species in MBW 
### NO DISTANCE SAMPLING
### and for all individuals detected - can back track + edit for heard vs seen

# switch to output folder

library(plyr)
library(unmarked)
library(ggplot2)
library(lubridate)
library(MuMIn)
library(plotrix)

setwd("C:/Users/akouw/OneDrive/Desktop/R/MBW2023/Output")

#sp.list<-c("bith","bcch","blpw","boch","fosp","wiwr", "resq", "ybfl", "heth", "swth", "wtsp")

#sp.list <-c("heth", "swth", "wtsp")


### 2023 data only

#sp.fn<-function(sp){
  
### to run comparison b/w 2012-2015 and 2016, change sp below to each species and run code
### until you get to models, then run only the model that has the highest AIC for each species
### then head to bottom of this page to run the abundance codes for each species


####2023: Run each species through code separately. Take best model and enter into code at bottom to get abundances.

sp<-"ybfl"

### SP PRESENT
pointID <- na.omit(pointID)

SP.PRES<-subset(all.sp, spcd==sp)

SP.PRES<-ddply(SP.PRES,.(pointID),function(df)
		data.frame(
		count.1 = sum(df$count.1),
		count.2 = sum(df$count.2),
		count.3 = sum(df$count.3),
		count.4 = sum(df$count.4)))

# In 2013 script have sight/sound/total count data
#	missing from 2012 data, so can just do sp.pres instead of
#	multiple divisions

sp.pres<-SP.PRES
names(sp.pres)<-c("pointID","t1","t2","t3","t4")

### SP ABSENT

sp.abs1<-subset(stop.dat,survey=="y",select=c(1:3))
names(sp.abs1)<-c("pointID","lat","long")

sp.abs2<-pointID[,c(1,9)]

sp.abs<-merge(sp.abs1,sp.abs2,by="pointID",all.x=T,all.y=T)

#sp.abs <- na.omit(sp.abs)

a<-as.data.frame(unique(sp.abs$period))
names(a)<-c("period")
a$t1<-0
a$t2<-0
a$t3<-ifelse(a$period=="am",0,NA)
a$t4<-ifelse(a$period=="am",0,NA)

all.abs<-merge(sp.abs,a,by=c("period"))

sp.abs<-all.abs[,c(2,5:8)]

rm(a)

## remove points from sp.abs where spp is actually present
sp.abs<-subset(sp.abs,!pointID %in% sp.pres$pointID)

sp.dat<-rbind(sp.pres,sp.abs)

## sort data to be sure it's in the same order for the 
## sp matrix + the Covs matricies

sp.dat<-sp.dat[with(sp.dat,order(pointID)),]

### Will have to add another step to split out by year in the future
### but for now am only 

### Should be able to set up UMF for all variables here ... 

### ALL DATA 10min

# 2023 data only:
sp.dat1 <- merge(sp.dat, pointID)
#sp.dat1 <- na.omit(sp.dat1)
sp.mat1<-as.matrix(sp.dat1[,c(2,3)])

# site covs variables
covs.1a<-merge(all.abs,pointID,by=c("pointID","period"))

covs.1a<-covs.1a[,c(1,3,4,2,9,11:19)]
covs.1a$routeID<-as.numeric(substr(covs.1a$pointID,6,7))
covs.1a$start<-covs.1a$time$hour+covs.1a$time$min/60



## 2023 data only 
# b/c only one site visit need to use pcount

sp.umf1<-unmarkedFramePCount(y=sp.mat1,
	siteCovs=data.frame(ptID=covs.1a$pointID,route=as.factor(covs.1a$routeID),
		lat=covs.1a$lat,long=covs.1a$long,
		year=covs.1a$year, obs=covs.1a$obs, jday=covs.1a$jday,
		start=covs.1a$start,time=covs.1a$period))


## this is where you stop for calculating abundances mentioned above and below

## null model
m0<-pcount(~1~1,data=sp.umf1,K=100,se=T)

## abundance models
m11<-pcount(~1 ~(start-1), data=sp.umf1,K=100,se=T)
m12<-pcount(~1 ~(start-1)/jday, data=sp.umf1,K=100,se=T)


## detection models
m13<-pcount(~(start-1) ~1, data=sp.umf1,K=100,se=T)
m14<-pcount(~jday ~1, data=sp.umf1,K=100,se=T)

## additive models
m15<-pcount(~1 ~(start-1+jday),data=sp.umf1,K=100,se=T)
m16<-pcount(~(start-1+jday) ~1, data=sp.umf1,K=100,se=T)
m17<-pcount(~(start-1+jday) ~(start-1+jday), data=sp.umf1,K=100,se=T)

## route specific models - likely crashy but can run or ## out if needed
#m11<-pcount(~1 ~(prov-1)/(route-1), data=sp.umf1,K=20,se=F)
#m12<-pcount(~(prov-1)/(route-1) ~1, data=sp.umf1,K=20,se=F)

## could just as easily use period (aka time) as a categorical variable but I prefer continuous

## get fit list + model covariates + model selection table
#a.fit<-fitList(m0,m13,m14)
#a.fit<-fitList(m0,m11,m12,m13,m14)
a.fit<-fitList(m0,m13,m14,m15,m16,m17)
#a.fit <-fitList(m0,m11,m12,m13,m14,m15,m16,m17)
ms.a<-modSel(a.fit) # model selection table

# Export model coefficiants
ce.a<-coef(ms.a) # model coefficients
ce.a$spcd<-sp
ce.a$mod<-row.names(ce.a)
write.table(ce.a,"modCoef.csv",sep=",", row.names=F,col.names=F,append=T)

# Export model SEs
se.a<-SE(ms.a) # model standard errors
se.a$spcd<-sp
se.a$mod<-row.names(se.a)
write.table(se.a,"modSE.csv",sep=",", row.names=F,col.names=F, append=T)

# Export AIC table
aic.a<-ms.a@Full[,c(1:2,14,16,17,18,20)] # AIC table
#aic.a$spcd<-sp
write.table(aic.a,"modAIC.csv",sep=",", row.names=F, col.names=F,append=T)


# Use model averaging to get population estimates (not the way I did it before)
# set up fake data (which is esentially the same as the original data)

x<-covs.1a[,c("pointID","routeID","lat","long", "day","month","year","jday","start","start_hour","start_min","period")]

x1<-ddply(x,.(routeID,period),function(df)min(df$start))
names(x1)[3]<-"start.min"

fake.dat<-merge(x,x1,by=c("routeID","period"))
fake.dat$d.start<-fake.dat$start-fake.dat$start.min
fake.dat$spcd<-sp
fake.dat$route<-as.factor(fake.dat$routeID)

rm(x1,x)

# get model predicted averages
# abundance
p1<-predict(a.fit,type="state",newdata=fake.dat,appendData=T)

# detection
p2<-predict(a.fit,type="det",newdata=fake.dat,appendData=T)

## export model predictions
write.table(p1,"modABU.csv",sep=",",row.names=F,col.names=F,append=T)
write.table(p2,"modDET.csv",sep=",",row.names=F,col.names=F,append=T)

#}
#lapply(sp.list,sp.fn)



#### *** Do the abundance model (code further below) for each species directly after you finish the above code for that species! *** 

##############################
##############################
##### Only do this part for first species through. Adds column labes in excel file.
# reimport data files

# AIC
aic.a<-read.csv("modAIC.csv",header=F)
colnames(aic.a)<-c("model","formula","nPars","AIC","delta","AICwt", "cumltvWt")
write.table(aic.a,"modAIC.csv",sep=",", row.names=F)

# Model coefficents
ce.a<-read.csv("modCoef.csv",header=F)
colnames(ce.a)<-c("lam(Int)","p(Int)","p(start)","p(jday)", "spcd","mod")
write.table(ce.a,"modCoef.csv",sep=",",row.names=F)

## Model SEs
## removed model SE's b/c otherwise the route specific models won't run (i.e. n==0)
## might work in subsequent years with more data ... but not sure

#se.a<-read.csv("modSE.csv",header=F)
#colnames(se.a)<-c("SElam(Int)","SElam(provNB)","SElam(provNB:jday)",
#	"SElam(provNB:start)","SElam(provNS)","SElam(provNS:jday)",
#	"SElam(provNS:start)","SEp(Int)","SEp(provNB)",       
#	"SEp(provNB:jday)","SEp(provNB:start)","SEp(provNS)",       
#	"SEp(provNS:jday)","SEp(provNS:start)","spcd","mod")
#write.table(se.a,"modSE.csv",sep=",", row.names=F,col.names=F, append=T)

#####################################
#####################################
### Only do this part for first species through. Adds column labes in excel file.
p1<-read.csv("modABU.csv",header=F)
colnames(p1)<-c("Predicted","SE","lower","upper","routeID", "period","pointID","lat","long",
	"day","month","year","jday","start", "start_hour","start_min","start.min","d.start","spcd","route")
write.table(p1,"modABU.csv",sep=",",row.names=F)

p2<-read.csv("modDET.csv",header=F)
colnames(p2)<-c("Predicted","SE","lower","upper","routeID", "period","pointID","lat","long",
	"day","month","year","jday","start",
	"start_hour","start_min","start.min","d.start","spcd","route")
write.table(p2,"modDET.csv",sep=",",row.names=F)



### Use code below to get abundance numbers for 2023.  Do each species separately above first and add in best model for each here.

## ABUNDANCE MODELS
pt<-covs.1a[,c(1,4,12)]

rebith<-ranef(m14)
bithab.est<-bup(rebith)
bithab.est<-cbind(pt,bithab.est)
bith.er <- t(std.error(bithab.est[,4]))
bith.mean <- mean(bithab.est[,4])
bithabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(bithabund) <- c("mean","se","year")
bithabund[1,2] <- t(bith.er)
bithabund[1,1] <- t(bith.mean)
bithabund[1,3] <- "2023"
write.table(bithabund,"abund2023.csv",sep=",",row.names=F)


rebcch<-ranef(m16)
bcchab.est<-bup(rebcch)
bcchab.est<-cbind(pt,bcchab.est)
bcch.er <- t(std.error(bcchab.est[,4]))
bcch.mean <- mean(bcchab.est[,4])
bcchabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(bcchabund) <- c("mean","se","year")
bcchabund[1,2] <- t(bcch.er)
bcchabund[1,1] <- t(bcch.mean)
bcchabund[1,3] <- "2023"
write.table(bcchabund,"abund2023.csv",sep=",",row.names=F)


reblpw<-ranef(m15)
blpwab.est<-bup(reblpw)
blpwab.est<-cbind(pt,blpwab.est)
blpw.er <- t(std.error(blpwab.est[,4]))
blpw.mean <- mean(blpwab.est[,4])
blpwabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(blpwabund) <- c("mean","se","year")
blpwabund[1,2] <- t(blpw.er)
blpwabund[1,1] <- t(blpw.mean)
blpwabund[1,3] <- "2023"
write.table(blpwabund,"abund2023.csv",sep=",",row.names=F)


reboch<-ranef(m16)
bochab.est<-bup(reboch)
bochab.est<-cbind(pt,bochab.est)
boch.er <- t(std.error(bochab.est[,4]))
boch.mean <- mean(bochab.est[,4])
bochabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(bochabund) <- c("mean","se","year")
bochabund[1,2] <- t(boch.er)
bochabund[1,1] <- t(boch.mean)
bochabund[1,3] <- "2023"
write.table(bochabund,"abund2023.csv",sep=",",row.names=F)


refosp<-ranef(m14)
fospab.est<-bup(refosp)
fospab.est<-cbind(pt,fospab.est)
fosp.er <- t(std.error(fospab.est[,4]))
fosp.mean <- mean(fospab.est[,4])
fospabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(fospabund) <- c("mean","se","year")
fospabund[1,2] <- t(fosp.er)
fospabund[1,1] <- t(fosp.mean)
fospabund[1,3] <- "2023"
write.table(fospabund,"abund2023.csv",sep=",",row.names=F)


reheth<-ranef(m14)
hethab.est<-bup(reheth)
hethab.est<-cbind(pt,hethab.est)
heth.er <- t(std.error(hethab.est[,4]))
heth.mean <- mean(hethab.est[,4])
hethabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(hethabund) <- c("mean","se","year")
hethabund[1,2] <- t(heth.er)
hethabund[1,1] <- t(heth.mean)
hethabund[1,3] <- "2023"
write.table(hethabund,"abund2023.csv",sep=",",row.names=F)


reswth<-ranef(m14)
swthab.est<-bup(reswth)
swthab.est<-cbind(pt,swthab.est)
swth.er <- t(std.error(swthab.est[,4]))
swth.mean <- mean(swthab.est[,4])
swthabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(swthabund) <- c("mean","se","year")
swthabund[1,2] <- t(swth.er)
swthabund[1,1] <- t(swth.mean)
swthabund[1,3] <- "2023"
write.table(swthabund,"abund2023.csv",sep=",",row.names=F)


rewiwr<-ranef(m15)
wiwrab.est<-bup(rewiwr)
wiwrab.est<-cbind(pt,wiwrab.est)
wiwr.er <- t(std.error(wiwrab.est[,4]))
wiwr.mean <- mean(wiwrab.est[,4])
wiwrabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(wiwrabund) <- c("mean","se","year")
wiwrabund[1,2] <- t(wiwr.er)
wiwrabund[1,1] <- t(wiwr.mean)
wiwrabund[1,3] <- "2023"
write.table(wiwrabund,"abund2023.csv",sep=",",row.names=F)


rewtsp<-ranef(m13)
wtspab.est<-bup(rewtsp)
wtspab.est<-cbind(pt,wtspab.est)
wtsp.er <- t(std.error(wtspab.est[,4]))
wtsp.mean <- mean(wtspab.est[,4])
wtspabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(wtspabund) <- c("mean","se","year")
wtspabund[1,2] <- t(wtsp.er)
wtspabund[1,1] <- t(wtsp.mean)
wtspabund[1,3] <- "2023"
write.table(wtspabund,"abund2023.csv",sep=",",row.names=F)


reybfl<-ranef(m16)
ybflab.est<-bup(reybfl)
ybflab.est<-cbind(pt,ybflab.est)
ybfl.er <- t(std.error(ybflab.est[,4]))
ybfl.mean <- mean(ybflab.est[,4])
ybflabund=data.frame(matrix(NA, nrow=1, ncol=3))
colnames(ybflabund) <- c("mean","se","year")
ybflabund[1,2] <- t(ybfl.er)
ybflabund[1,1] <- t(ybfl.mean)
ybflabund[1,3] <- "2023"
write.table(ybflabund,"abund2023.csv",sep=",",row.names=F)



spcd <- c("bith","bcch","blpw","boch","fosp","heth","swth","wiwr","wtsp","ybfl")
abund2023 <- rbind(bithabund, bcchabund, blpwabund, bochabund, fospabund, hethabund, swthabund, wiwrabund, wtspabund, ybflabund)
abund2023 <- cbind(spcd,abund2023)
write.table(abund2023,"allabund2023.csv",sep=",",row.names=F)

