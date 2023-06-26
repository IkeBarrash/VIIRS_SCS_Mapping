#Ike Barrash

#The Chinese Maritime Militia (CMM) are a paramilitary force used by the People's Republic of China
#further Chinese maritime interests through grey-zone actions--activities which do not cross thresholds
#that would result in a military response. CMM Ships are numerous, and difficult to track using 
#traditional ship tracking techniques, and traverse large swaths of oceans making manual detection 
#unfeasible.

#Detecting and countering grey-zone activities has garnered increasing attention in U.S. talks with
#regional countries. It is often not possible to share classified information the U.S. gathers on the
#region creating a market for open-source intelligence solutions.

#This code utilizes automated detection of ship lights at night from the South China Sea region drawn 
#from the School of Mines VIIRS database to determine how many ships are located around a list of South
#China Sea features in the contested Spratly Island and Paracel Island groups. This code creates two
#sets of output. First daily observations are plotted on a video map. Next, ship counts within a certain
#distance of each feature over time are graphed, as well as overall ship counts.

#This program allows quick comparison of ship numbers around a large number of ship features over time.
#Using this data allows the user to quickly determine where grey-zone activities may have occurred, and
#project potential locations and dates to investigate further.


rm(list=ls())

#Libraries needed
library(tidyverse)
library(ggplot2)
library(lubridate)
library(data.table)
library(imputeTS)
library(maps)
library(mapdata)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(gganimate)
library(gifski)
library(readxl)
library(sp)
#Making a list of all VIIRS files by looking for "VBD" in the file name
filenames = list.files(path = "C:/Users/ikeba/Desktop/Work/CSIS/VIIRS", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
filenames = filenames[grep(pattern = "VBD",filenames)]

#Setting file path to folder with downloaded VIIRS daily CSV files ###Change for other users###
#files available at https://eogdata.mines.edu/products/vbd/#download
setwd("C:/Users/ikeba/Desktop/Work/CSIS/VIIRS")

#Reading detections of each day and appending to a master dataset
for(i in 1:length(filenames)){
 toadd = fread(input = filenames[i],header=T,select=c("id_Key","Lat_DNB","Lon_DNB","Date_Mscan","QF_Detect","Land_Mask"))
 toadd = toadd%>%mutate(id_Key=as.character(id_Key),
                        Lat_DNB = as.numeric(Lat_DNB),
                        Lon_DNB = as.numeric(Lon_DNB),
                        Date_Mscan = ymd(as.character(floor_date(ymd_hms(Date_Mscan), unit = "day"))), #Making date variable and removing time from date
                        QF_Detect = as.numeric(QF_Detect),
                        Land_Mask = as.numeric(Land_Mask))
  if(i == 1){
    Dat = toadd
  }
  else{
    Dat = rbind(Dat,toadd)
  }
}

#Removing NAs (blank observation at the end of each file is recorded as an NA)
summary(Dat)
Dat = na.omit(Dat)
Dat = Dat%>%filter(Lat_DNB<100)

#Loading HD world Data
world = ne_countries(scale = "large", returnclass = "sf")

#Filtering out all non-ship type detections
Dat = Dat%>%filter(QF_Detect == (1|2))

#######################Mapping observations over time##########################
#Each map runs from October 2021-June 2022 and produces a moving map of where 
#ships were seen each day

#SCS MAP
SCSdaymap = ggplot(data = world) +
  geom_sf(color = "black", fill = "darkgreen") +
  coord_sf(xlim = c(107, 122.5), ylim = c(4, 19.5), expand = FALSE) +
  geom_point(data = Dat, aes(y=Lat_DNB,x=Lon_DNB),alpha = .1, color = "red") + 
  labs(title = 'South China Sea - VIIRS Detection Date: {frame_time}', x = "Longitude", y = "Latitude")+ theme(text = element_text(size = 20),panel.background = element_rect(fill = "lightblue"))+
  transition_time(Date_Mscan)
animate(SCSdaymap, height = 1200, width =1200,fps = 3,nframes = length(unique(Dat$Date_Mscan)))
anim_save('plot_VIIRS_data_Oct-June_SCS.gif')

#Zoomed in Spratlys map
SCSdaymap = ggplot(data = world) +
  geom_sf(color = "black", fill = "darkgreen") +
  coord_sf(xlim = c(110, 117), ylim = c(6.5, 14), expand = FALSE) +
  geom_point(data = Dat, aes(y=Lat_DNB,x=Lon_DNB),alpha = .1, color = "red") + 
  labs(title = 'Spratly Islands - VIIRS Detection Date: {frame_time}', x = "Longitude", y = "Latitude")+ theme(text = element_text(size = 20),panel.background = element_rect(fill = "lightblue"))+
  transition_time(Date_Mscan)
animate(SCSdaymap, height = 1200, width =1200,fps = 3,nframes = length(unique(Dat$Date_Mscan)))
anim_save('plot_VIIRS_data_Oct-June_Spratly.gif')

############Analyzing Observations Around South China Sea Features##############
#Defining distance within how many miles a detection can be to be associated with a feature
#(Can be changed based on preference of analyst)
Dist = 10
#Reading in AMTI Island database
Feat_Data = read_excel("C:/Users/ikeba/Desktop/Work/CSIS/SCS Occupied Features.xlsx")%>%select("name","latitude","longitude","occupier")

#Creating a data table with an observation for each feature on each day
Comp_Feat_Data = data.table("Feature" = sort(rep(unique(Feat_Data$name),length(min(Dat$Date_Mscan):max(Dat$Date_Mscan)*length(unique(Feat_Data$name))))),
                            "Date"= rep(min(Dat$Date_Mscan):max(Dat$Date_Mscan),length(unique(Feat_Data$name))))

#Merging with AMTI feature database
Comp_Feat_Data = Comp_Feat_Data%>%left_join(Feat_Data, by=c("Feature"="name"))

#Changing longitude and latitude to radians
Comp_Feat_Data = Comp_Feat_Data%>%mutate(lat1 = latitude/180*pi,long1 = longitude/180*pi)
Dat = Dat%>%mutate(lat1 = Lat_DNB/180*pi,long1 = Lon_DNB/180*pi)

###Checking every observation in Dat to find distance to each feature

#Defining function that checks if great circle distance between two points is below a distance D
GCDistCheck = function(lon1,lat1,lon2,lat2,D){
  if_else(D<3963.0*acos((sin(lat1)*sin(lat2))+cos(lat1)*cos(lat2)*cos(lon2-lon1)),0,1)
}
#Initializing count of ships to 0 for each day by each feature
Comp_Feat_Data = Comp_Feat_Data%>%mutate(Ships = 0)

#Loop counts number of ship detections within D distance of each feature on each day
for(i in 1:nrow(Comp_Feat_Data)){
  temp = Dat%>%filter(Date_Mscan == Comp_Feat_Data$Date[i])
  temp = temp%>%mutate(prox = GCDistCheck(long1,lat1,Comp_Feat_Data$long1[i],Comp_Feat_Data$lat1[i],Dist))
  Comp_Feat_Data$Ships[i] = sum(temp$prox)
}

#Recoding Occupier variable for legibility
Comp_Feat_Data = Comp_Feat_Data%>%mutate(occupier = factor(occupier))
Comp_Feat_Data$occupier = Comp_Feat_Data$occupier%>%fct_recode("China"="cn","Malaysia"="mys","Philippines"="php","Taiwan"="twn","Vietnam"="vn")

###############Plotting total number of ships over time colored by feature#############################
#Single stacked area chart
ggplot(Comp_Feat_Data) + geom_area(aes(x=as_date(Date),y=Ships,fill=Feature)) + 
  labs(x="Time", y="Number of Ships") + 
  scale_x_date(date_labels = "%b",date_breaks = "1 month") +
  scale_y_continuous(breaks = seq(0, 250, by = 20)) + theme_bw() + theme(plot.background = element_rect(fill = 'gray90'))

#Line charts faceted by Claimant
ggplot(Comp_Feat_Data) + geom_line(aes(x=as_date(Date),y=Ships,color=Feature)) + 
  labs(x="Time", y="Number of Ships") + 
  scale_x_date(date_labels = "%b",date_breaks = "1 month") +
  scale_y_continuous(breaks = seq(0, 100, by = 4)) + theme_bw() + facet_wrap(vars(occupier)) + theme(plot.background = element_rect(fill = 'gray90'))

#Future work:
#
#
#Write code that throws out observations that aren't within a certain distance of any feature,
#Create separate data.table which only counts ships for the island they are closest too as
#opposed to all withing D.



