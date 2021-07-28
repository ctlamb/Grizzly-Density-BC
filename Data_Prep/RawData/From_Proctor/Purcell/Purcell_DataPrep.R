#################################
## PREP PROCTOR CAPTURE DATA FOR USE IN DENSITY META-ANALYSIS
## PURCELLS
#################################

##Load Spatial packages
require(sp)
require(rgdal)
require(maps)
require(raster)
require(rgeos)
require(maptools)
require(RgoogleMaps)
require(ggmap)
require(plyr)

##################################################################
####LOAD CAPTURE DATA
##################################################################

#set working directory
setwd("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_Proctor")


cap <- read.csv("Purcell/1Purcell Capture data for JB_edited_CL.csv", header=TRUE, sep=",")
site <- read.csv("1SS Pur DNA sites all_edited_CL.csv", header=TRUE, sep=",") 

site <-rename(site ,c("YEAR"="Year"))
##################################################################
##### PLOT MAP ALL OF ALL SITE DATA
##################################################################

##Cull any missing spatial data
site<-subset(site,!is.na(EAST))

####Make Data Spatial
coordinates(site) <- c("EAST","NORTH")
proj4string(site)<-"+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"



##Cull any missing data

site$Easting.DC<-coordinates(site)[,1]
site$Northing.DC<-coordinates(site)[,2]

##transform to lat long
site<-spTransform(site,CRS("+init=epsg:4326"))
site<-as.data.frame(site)

###Plot to make sure all is good

basemap <- get_map(location = c(lon = -116.5, lat = 49.7), zoom = 8)
map <- ggmap(basemap)+
  geom_point(aes(x = EAST,
                 y = NORTH,
                 fill=PROJECT),
             colour = "black",
             pch = 21,
             data=site,size=5)

map

##################################################################
##### SELECT ONLY PROJECT OF INTEREST FROM SITE DATA AND LOOP

##FROM PROCTOR:

#####Keep in mind, a few of these surveys were complete grids.
#SS 2005
#Hwy 3 East 2004
#Hwy 3 West 2005
#Jumbo 1998

##The rest were done survey like, that is they had 4 2 week collections session, but were just put in the best of sites for my fragmentation PhD work.
##I did it that way, so one day I could hobble a pop est out of the data.
#SP 98
#SP 2001
#PWC 2002

##And 2 where i just did a few sites to catch a few extra bears. Don't use these.
#SP 99
#NP 99

##################################################################

##Rename Studies to More Descriptive Names
cap$Project<-as.character(cap$Project)
cap[cap$Project=="Hwy3 west", "Project"]<-"HWY3_2005"
cap[cap$Project=="Hwy3 east", "Project"]<-"HWY3_2004"
cap[cap$Project=="PWC", "Project"]<-"PWC_2002"
cap[cap$Project=="SP01", "Project"]<-"S_Purcell_2001"
cap[cap$Project=="SP98", "Project"]<-"S_Purcell_1998"
cap[cap$Project=="Jumbo", "Project"]<-"Jumbo_1998"


##Rename so site file aligns
site$PROJECT<-as.character(site$PROJECT)
site[site$PROJECT=="Hwy3 2005", "PROJECT"]<-"HWY3_2005"
site[site$PROJECT=="Hwy 3 2004", "PROJECT"]<-"HWY3_2004"
site[site$PROJECT=="PWC 2002", "PROJECT"]<-"PWC_2002"
site[site$PROJECT=="Purcell 2001", "PROJECT"]<-"S_Purcell_2001"
site[site$PROJECT=="SP 1998", "PROJECT"]<-"S_Purcell_1998"
site[site$PROJECT=="Jumbo 1998", "PROJECT"]<-"Jumbo_1998"


##Drop studies from site data that are not going to join
site<-subset(site, !PROJECT %in% c("SP 1999", "N Purcell 1999","S Selkirk 2005"))



###Plot to make sure all is good

basemap <- get_map(location = c(lon = -116.5, lat = 49.7), zoom = 8)
map <- ggmap(basemap)+
  geom_point(aes(x = EAST,
                 y = NORTH,
                 fill=PROJECT),
             colour = "black",
             pch = 21,
             data=site,size=5)

map



##################################################################
##### SET UP LOOP
##################################################################

Proj<-as.character(unique(cap$Project))

for(i in 1:length(Proj)){
  
  
  ##################################################################
  ##### PREP CAPTURE DATA
  ##################################################################
  
  ##Subset out project of interest
  ##retain useful columns
  cap_join <- cap[cap$Project==Proj[i],c("Year", "Perm.ID", "Sex", "Site", "Session", "Site.name", "Easting", "Northing")]
  

  ###Change sessions from letters to #'s, A's are equlivelent to sess 1, B=2, C=3, D=4
  cap_join$Session<-as.character(cap_join$Session)
  cap_join$Session <-ifelse( cap_join$Session == "A", as.character(1), cap_join$Session )
  cap_join$Session <-ifelse( cap_join$Session == "B", as.character(2), cap_join$Session )
  cap_join$Session <-ifelse( cap_join$Session == "C", as.character(3), cap_join$Session )
  cap_join$Session <-ifelse( cap_join$Session == "D", as.character(4), cap_join$Session )
  
  

##################################################################
##### PREP SITE DATA
##################################################################
proj_site<- site[ site$PROJECT==Proj[i],]
  
  
##Ensure no duplicates
##Pull out unique records, based on certain columns while retaining all columns
proj_site <-  proj_site[!duplicated( proj_site[,c( "Year", "SITE_NAME", "Easting.DC", "Northing.DC" )]),]

###Plot to make sure all is good

basemap <- get_map(location = c(lon = -116.5, lat = 49.7), zoom = 8)
map <- ggmap(basemap)+
  geom_point(aes(x = EAST,
                 y = NORTH,
                 fill=PROJECT),
             colour = "black",
             pch = 21,
             data=proj_site,size=5)

print(map)

ggsave(paste("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/Maps/",Proj[i],".pdf",sep=""))

##retain useful columns
proj_site <- proj_site[,c("Year","SITE_NAME", "Easting.DC", "Northing.DC")]





if( max(as.numeric(cap_join$Session)) != 4) stop("Sessions != 4")

###These sites were deployed for four sessions and not moved
##duplicate sites 4x, with session # 1-4 each time
site_join<-data.frame(rbind(proj_site,
                            proj_site,
                            proj_site,
                            proj_site),
                      Session=rep(c("1","2","3","4"),each=nrow(proj_site)))








##################################################################
##### MERGE DATA
##################################################################

##merge site and capture data
FINAL<-merge(site_join, cap_join, by.x=c("Year","SITE_NAME", "Session", "Easting.DC","Northing.DC"),  by.y=c("Year","Site.name", "Session", "Easting", "Northing"), all=TRUE)

###checks
if(length(unique(FINAL$SITE_NAME))!=length(unique(site_join$SITE_NAME))) warning(paste(Proj[i], "Not same # of Sites, Match names",sep=" "))

if(length(unique(FINAL$SITE_NAME))!=length(unique(site_join$SITE_NAME))) { 
  print(unique(FINAL$SITE_NAME)[!unique(FINAL$SITE_NAME)%in%unique(site_join$SITE_NAME)])
  print(unique(site_join$SITE_NAME[order(site_join$SITE_NAME)]))
  }

if(length(unique(FINAL[!is.na(FINAL$Perm.ID),]$Perm.ID))!=length(unique(cap_join$Perm.ID))) warning(paste(Proj[i], "Not same # of Individuals, Good luck",sep=" "))

if(nrow(FINAL[!is.na(FINAL$Perm.ID),])!=nrow(cap_join)) warning(paste(Proj[i], "Not same # of Individuals #2, Good luck",sep=" "))



##Add Trap Type
FINAL$Trap.Type<-"BS"

##Add Project
FINAL$Project<-as.character(Proj[i])

##Add Zone
FINAL$Zone<-11

##Rename
FINAL<-rename(FINAL,c("Perm.ID"="Individual"))

##Create unique Site ID for each site
FINAL$SITE_NAME<-paste(FINAL$SITE_NAME,FINAL$Easting.DC,FINAL$Northing.DC,sep="_")

##Reorder
FINAL <- FINAL[, c("Project","Year", "SITE_NAME", "Trap.Type",  "Session",     "Easting.DC",  "Northing.DC", "Zone", "Individual"  , "Sex"   )]
FINAL<-FINAL[with(FINAL, order(SITE_NAME, Session)), ]

##EXPORT
write.csv(FINAL, paste("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/ProjectFiles/Single_Year/",Proj[i], ".csv", sep=""), row.names = FALSE, na="")

}

