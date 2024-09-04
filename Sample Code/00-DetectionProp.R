# Calculate detection probabilities across years for CWMP frog trends paper
# writen by Doug Tozer 

library(ggplot2)
library(png)
library(dplyr)
library(tidyr)
library(unmarked)

# Bring in raw data

data <- read.csv("Data/CWMP_anuran_obs.csv")
summary(data)
str(data)
names(data)
head(data)
nrow(data)

# Create a factor version of year for later use (year2)
# And make some new variables for later to match Unmarked verbage

data <- data %>%
	mutate(year2=as.factor(year)) %>%
	mutate(M=paste(site_id,point_id,year,sep="-")) %>%
	mutate(J=ifelse(sample==1,"visit1",ifelse(sample==2,"visit2","visit3")))

# Create speciesMaster

speciesMaster <- data %>%
	filter(taxa_code!="(none)") %>%
	group_by(M,year2,J,taxa_code) %>%
	summarize() %>%
	mutate(occurrence=1)

# Create newdata for use later

newdata <- data.frame(year2=seq(2011,2023,1))
newdata$year2 <- as.factor(newdata$year2)

#################################
AMTO
#################################

amtoTemp <- speciesMaster %>%
	filter(taxa_code=="AMTO")

amtoY <- amtoTemp %>%
	spread(J,occurrence) %>%
	replace(is.na(.), 0) %>%
	ungroup() %>%
	arrange(M) %>%
	select(visit1,visit2,visit3)

amtoSiteCovs <- amtoTemp %>%
	group_by(M,year2) %>%
	summarize() %>%
	ungroup() %>%
	arrange(M) %>%
	select(year2)
	
amtoOccuData <- unmarkedFrameOccu(amtoY,siteCovs=amtoSiteCovs)

summary(amtoOccuData)	

system.time(amtoOccuModel <- occu(~year2 ~1, data=amtoOccuData))

amtoOccuPreds <- data.frame(species="amto",predict(amtoOccuModel, type="det", newdata=newdata, appendData=TRUE))

#################################
BULL
#################################

bullTemp <- speciesMaster %>%
	filter(taxa_code=="BULL")

bullY <- bullTemp %>%
	spread(J,occurrence) %>%
	replace(is.na(.), 0) %>%
	ungroup() %>%
	arrange(M) %>%
	select(visit1,visit2,visit3)

bullSiteCovs <- bullTemp %>%
	group_by(M,year2) %>%
	summarize() %>%
	ungroup() %>%
	arrange(M) %>%
	select(year2)
	
bullOccuData <- unmarkedFrameOccu(bullY,siteCovs=bullSiteCovs)

summary(bullOccuData)	

system.time(bullOccuModel <- occu(~year2 ~1, data=bullOccuData))

bullOccuPreds <- data.frame(species="bull",predict(bullOccuModel, type="det", newdata=newdata, appendData=TRUE))

#################################
CHFR
#################################

chfrTemp <- speciesMaster %>%
	filter(taxa_code=="CHFR")

chfrY <- chfrTemp %>%
	spread(J,occurrence) %>%
	replace(is.na(.), 0) %>%
	ungroup() %>%
	arrange(M) %>%
	select(visit1,visit2,visit3)

chfrSiteCovs <- chfrTemp %>%
	group_by(M,year2) %>%
	summarize() %>%
	ungroup() %>%
	arrange(M) %>%
	select(year2)
	
chfrOccuData <- unmarkedFrameOccu(chfrY,siteCovs=chfrSiteCovs)

summary(chfrOccuData)	

system.time(chfrOccuModel <- occu(~year2 ~1, data=chfrOccuData))

chfrOccuPreds <- data.frame(species="chfr",predict(chfrOccuModel, type="det", newdata=newdata, appendData=TRUE))

#################################
GRTR
#################################

grtrTemp <- speciesMaster %>%
	filter(taxa_code=="GRTR")

grtrY <- grtrTemp %>%
	spread(J,occurrence) %>%
	replace(is.na(.), 0) %>%
	ungroup() %>%
	arrange(M) %>%
	select(visit1,visit2,visit3)

grtrSiteCovs <- grtrTemp %>%
	group_by(M,year2) %>%
	summarize() %>%
	ungroup() %>%
	arrange(M) %>%
	select(year2)
	
grtrOccuData <- unmarkedFrameOccu(grtrY,siteCovs=grtrSiteCovs)

summary(grtrOccuData)	

system.time(grtrOccuModel <- occu(~year2 ~1, data=grtrOccuData))

grtrOccuPreds <- data.frame(species="grtr",predict(grtrOccuModel, type="det", newdata=newdata, appendData=TRUE))

#################################
GRFR
#################################

grfrTemp <- speciesMaster %>%
	filter(taxa_code=="GRFR")

grfrY <- grfrTemp %>%
	spread(J,occurrence) %>%
	replace(is.na(.), 0) %>%
	ungroup() %>%
	arrange(M) %>%
	select(visit1,visit2,visit3)

grfrSiteCovs <- grfrTemp %>%
	group_by(M,year2) %>%
	summarize() %>%
	ungroup() %>%
	arrange(M) %>%
	select(year2)
	
grfrOccuData <- unmarkedFrameOccu(grfrY,siteCovs=grfrSiteCovs)

summary(grfrOccuData)	

system.time(grfrOccuModel <- occu(~year2 ~1, data=grfrOccuData))

grfrOccuPreds <- data.frame(species="grfr",predict(grfrOccuModel, type="det", newdata=newdata, appendData=TRUE))

#################################
NLFR
#################################

nlfrTemp <- speciesMaster %>%
	filter(taxa_code=="NLFR")

nlfrY <- nlfrTemp %>%
	spread(J,occurrence) %>%
	replace(is.na(.), 0) %>%
	ungroup() %>%
	arrange(M) %>%
	select(visit1,visit2,visit3)

nlfrSiteCovs <- nlfrTemp %>%
	group_by(M,year2) %>%
	summarize() %>%
	ungroup() %>%
	arrange(M) %>%
	select(year2)
	
nlfrOccuData <- unmarkedFrameOccu(nlfrY,siteCovs=nlfrSiteCovs)

summary(nlfrOccuData)	

system.time(nlfrOccuModel <- occu(~year2 ~1, data=nlfrOccuData))

nlfrOccuPreds <- data.frame(species="nlfr",predict(nlfrOccuModel, type="det", newdata=newdata, appendData=TRUE))

#################################
SPPE
#################################

sppeTemp <- speciesMaster %>%
	filter(taxa_code=="SPPE")

sppeY <- sppeTemp %>%
	spread(J,occurrence) %>%
	replace(is.na(.), 0) %>%
	ungroup() %>%
	arrange(M) %>%
	select(visit1,visit2,visit3)

sppeSiteCovs <- sppeTemp %>%
	group_by(M,year2) %>%
	summarize() %>%
	ungroup() %>%
	arrange(M) %>%
	select(year2)
	
sppeOccuData <- unmarkedFrameOccu(sppeY,siteCovs=sppeSiteCovs)

summary(sppeOccuData)	

system.time(sppeOccuModel <- occu(~year2 ~1, data=sppeOccuData))

sppeOccuPreds <- data.frame(species="sppe",predict(sppeOccuModel, type="det", newdata=newdata, appendData=TRUE))

#################################
WOFR
#################################

wofrTemp <- speciesMaster %>%
	filter(taxa_code=="WOFR")

wofrY <- wofrTemp %>%
	spread(J,occurrence) %>%
	replace(is.na(.), 0) %>%
	ungroup() %>%
	arrange(M) %>%
	select(visit1,visit2,visit3)

wofrSiteCovs <- wofrTemp %>%
	group_by(M,year2) %>%
	summarize() %>%
	ungroup() %>%
	arrange(M) %>%
	select(year2)
	
wofrOccuData <- unmarkedFrameOccu(wofrY,siteCovs=wofrSiteCovs)

summary(wofrOccuData)	

system.time(wofrOccuModel <- occu(~year2 ~1, data=wofrOccuData))

wofrOccuPreds <- data.frame(species="wofr",predict(wofrOccuModel, type="det", newdata=newdata, appendData=TRUE))



##################################################
# Graph it!                                      #
##################################################

predData <- rbind(amtoOccuPreds,bullOccuPreds,chfrOccuPreds,grtrOccuPreds,grfrOccuPreds,nlfrOccuPreds,sppeOccuPreds,wofrOccuPreds)

relabel <- data.frame(species=c("amto","bull","chfr","grtr","grfr","nlfr","sppe","wofr"),
		species2=c("American Toad","Bullfrog","Chorus Frog","Common Gallinule","Common Grackle","Common Yellowthroat","Forster's Tern","Least Bittern", "Marsh Wren","Mute Swan","Pied-billed Grebe","Red-winged Blackbird","Sandhill Crane","Sedge Wren","Sora","Swamp Sparrow","Virginia Rail","Wilson's Snipe"))

predData <- predData %>%
left_join(relabel)

detPlot <- ggplot(data=predData,aes(x=year2,y=Predicted))+
	geom_point(shape=21)+
	geom_errorbar(aes(ymin = lower, ymax = upper))+
	geom_smooth(aes(x=as.numeric(year2),y=Predicted),se=F)+
	labs(y="Detection probability")+
	ylim(0,1)+
	scale_x_discrete(breaks=c(2011,2015,2019))+
#	annotate(
	facet_wrap(~species2,ncol=3)+
	theme_bw() + theme(panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	plot.title = element_text(hjust = 0.5),axis.title.y=element_text(colour="black",size=12),
	axis.text=element_text(colour="black",size=10),legend.position = "top",
	axis.title.x=element_blank(),strip.text.x = element_text(size = 9.5),
	legend.text=element_text(size=12),legend.title=element_blank())

ggsave(file="1. det prob cwmp trends paper.png",plot=detPlot,width=8,height=10,units="in",dpi=800,
		path="C:/Users/dtozer/Documents/EPA Coastal/Manuscripts/Trends/Detection probability")




