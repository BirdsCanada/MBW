## Mountain Birdwatch Data Import #############################################


# Set up work space + load appropriate packages etc

#setwd(choose.dir())
setwd("C:/Users/akouw/OneDrive/Desktop/R/MBW2023") 

library(plyr)
library(unmarked)
library(ggplot2)
library(lubridate)
library(MuMIn)

### Data import - site covariates for MBW data #################################

# Survey data 
station<-read.csv("MBW2023NB_station.csv",header=T, sep=",")

### Stop ID data
stop <-read.csv("MBW2023NB_route.csv",header=T, sep=",")

## remove stops that weren't surveyed
stop.23<-subset(stop, surveyed.2023=="y")
stop.23<-stop.23[,c(1:4)]
stop.23$year<-2023

stop.dat<-stop[,c(1,3,4,29)]
names(stop.dat)[4]<-"survey.23"

### ALL MBW BIRD DATA ################################################
#
#	Analysis for all species in MBW (bith + fosp etc)
#	Data genates 3 species data tables
#		- visual detections
#		- auditory detections
#		- total detections

### 2023 DATA
bird.23<-read.csv("MBW2023NB_birds.csv",header=T, sep=",")


# JULIAN DATE for time sequence analysis
bird.23$date<-paste(bird.23$day,bird.23$month, bird.23$year,sep="/")
bird.23$date<-as.Date(bird.23$date, format='%d/%m/%Y')
bird.23$jday<-as.numeric(format(bird.23$date, format="%j"))


bird <- bird.23

table(bird$species)
table(bird$total_count)


#### Reformatting things ...


x1<-subset(bird, period==1)
x1<-x1[,c(1,17,7,15,9:11,8)]
names(x1)<-c("pointID","jday","year","period","spcd","min.d","max.d",
             "count.1")

x2<-subset(bird, period==2)
x2<-x2[,c(1,17,7,15,9:11,8)]
names(x2)<-c("pointID","jday","year","period","spcd","min.d","max.d",
             "count.2")

x3<-subset(bird, period==3)
x3<-x3[,c(1,17,7,15,9:11,8)]
names(x3)<-c("pointID","jday","year","period","spcd","min.d","max.d",
             "count.3")

x4<-subset(bird, period==4)
x4<-x4[,c(1,17,7,15,9:11,8)]
names(x4)<-c("pointID","jday","year","period","spcd","min.d","max.d",
             "count.4")


# match to the master data table ...

names(stop.dat)[1]<-"pointID"
names(stop.dat)[4]<-"survey"

b1<-merge(x1,x2,by=c("pointID","jday","year","period",
	"spcd","min.d","max.d"),all.x=T,all.y=T)

b2<-merge(b1,x3,by=c("pointID","jday","year","period",
	"spcd","min.d","max.d"),all.x=T,all.y=T)

sp.ALL<-merge(b2,x4,by=c("pointID","jday","year","period",
	"spcd","min.d","max.d"),all.x=T,all.y=T)

sp.ALL<-merge(stop.dat,sp.ALL,by=c("pointID"),all.x=T)

rm(b1,b2,x1,x2,x3,x4)

### Now add 0s to all species data
#  	- For AM + PM counts: NA in period 1 or 2 == 0
#	- For AM counts: NA in period 3 or 4 == 0
# 	- For PM counts: NA in period 3 of 4 == NA 

# sites not surveyed in a year 
x1<-subset(sp.ALL,survey=="n")

# all survyed sites in a year
x2<-subset(sp.ALL,survey=="y")
x2[,c(11:12)][is.na(x2[,c(11:12)])]<-0

# morning surveys (replace na with 0)
x3<-subset(x2,period=="am")
x3[,c(13:14)][is.na(x3[,c(13:14)])]<-0

# evening surveys
x4<-subset(x2,period=="pm")

# reassemble into mega data set

sp.all<-rbind(x1,x3,x4)

# sort it so it lines up "properly"
sp.all<-sp.all[with(sp.all,order(pointID)),]


### Dealing with distance sampling data now - making 2 data sets..

sp.DS<-sp.all

### No distance sampling (i.e. lump all distance categories)

all.sp<-ddply(sp.all,.(pointID,jday,year,spcd,period),function(df)
	data.frame(
		count.1 = sum(df$count.1),
		count.2 = sum(df$count.2),
		count.3 = sum(df$count.3),
		count.4 = sum(df$count.4)))

### SET UP SURVEY INFO (aka siteCovs + obsCovs df)
### basic pointID data 


y<-station[,c(1,3,6:9,21,22,23)]
names(y)[1]<-"pointID"
names(y)[2]<-"obs"
names(y)[3]<-"prov"

### GET a "master time data" column
y$chrono<-with(y,ymd_hm(paste(year,month,day,start_hour,start_min)))

### GET Julian day
y$jday<-as.numeric(format(y$chrono, format="%j"))

### GET start time
y$start<-format(y$chrono, format="%H:%M")
y$time<-strptime(y$start, format='%H:%M')

### POINT ID

pointID<-y
pointID <-na.omit(pointID)
#pointID <- pointID2

rm(y)






##Export for ACCDC format

bith.23<-read.csv("MBW2023NB_bith.csv",header=T, sep=",")
names(bith.23)[1]<-"pointID"
bith.23<-bith.23[,c(1:7)]


#bird.22.1<-bird.22[,c(1:11)]
#names(bird.22.1)[1]<-"pointID"
#bird.22.1<-subset(bird.22.1, species=="bith")
#bird.22.1<-bird.22.1[,c(1,11)]
#bird.22.2 <-subset(bird.22.1, pointID = bith.22)

pointID.1<-pointID[,c(1:5)]
names(pointID.1)[1]<-"pointID"

MBW23merge<- merge(bith.23, pointID.1,by=c("pointID"),all.x=T)
MBW23merge1<- merge(MBW23merge,stop.dat, by=c("pointID"),all.x=T)

write.csv(MBW23merge1, "MBW23_ACCDC.csv") 


### now should be able to switch scripts to do the sp loop

