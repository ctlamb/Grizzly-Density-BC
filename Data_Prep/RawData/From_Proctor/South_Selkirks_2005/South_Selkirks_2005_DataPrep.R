#################################
## PREP PROCTOR CAPTURE DATA FOR USE IN DENSITY META-ANALYSIS
## SOUTH SELKIRKS 2005
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

##################################################################
####LOAD CAPTURE DATA
##################################################################

#set working directory
setwd("/Users/clayton.lamb/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_Proctor")


cap <- read.csv("South_Selkirks_2005/1SS 2005 DNA results.csv", header=TRUE, sep=",")
sex <- read.csv("South_Selkirks_2005/1SS 2005 DNA results_SEX.csv", header=TRUE, sep=",")
site <- read.csv("1SS Pur DNA sites all_edited_CL.csv", header=TRUE, sep=",") 
bb <- read.csv("South_Selkirks_2005/w9845 S Selirks 2005 Results_CLEdits.csv", header=TRUE, sep=",") 

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
##### SELECT ONLY PROJECT OF INTEREST FROM SITE DATA
##################################################################

proj_site <- subset(site, PROJECT=="S Selkirk 2005")


###Plot to make sure all is good

basemap <- get_map(location = c(lon = -116.5, lat = 49.7), zoom = 8)
map <- ggmap(basemap)+
  geom_point(aes(x = EAST,
                 y = NORTH,
                 fill=PROJECT),
             colour = "black",
             pch = 21,
             data=proj_site,size=5)

map

ggsave("/Users/clayton.lamb/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/Maps/S_Selkirks_2005.pdf")

##Zoom in a bit
basemap <- get_map(location = c(lon = -117, lat = 49.35), zoom = 9)
map <- ggmap(basemap)+
  geom_point(aes(x = EAST,
                 y = NORTH,
                 fill=PROJECT),
             colour = "black",
             pch = 21,
             data=proj_site,size=5)

map



##################################################################
##### PREP SITE DATA
##################################################################

colnames(proj_site)

##retain useful columns

proj_site <- proj_site[,c("SITE_NAME", "Easting.DC", "Northing.DC")]

##Ensure no duplicates
##Pull out unique records, based on certain columns while retaining all columns
proj_site <- proj_site[!duplicated(proj_site[,c("SITE_NAME", "Easting.DC", "Northing.DC" )]),]

###These sites were deployed for four sessions and not moved
##duplicate sites 4x, with session # 1-4 each time
site_join<-data.frame(rbind(proj_site,
                            proj_site,
                            proj_site,
                            proj_site),
                      Session=rep(c("1","2","3","4"),each=nrow(proj_site)))




##################################################################
##### PREP CAPTURE DATA
##################################################################


##retain useful columns
cap_join <- cap[,c( "Individual", "Site.Sess", "Name", "East", "North")]

##seperate out session letter
cap_join$Site.Sess <- as.character(cap_join$Site.Sess)
cap_join$Session_let <- substr(cap_join$Site.Sess, start=nchar(cap_join$Site.Sess), stop=nchar(cap_join$Site.Sess)+1)

###Change sessions from letters to #'s, A's are equlivelent to sess 1, B=2, C=3, D=4
cap_join$Session <-NA
cap_join$Session <-ifelse( cap_join$Session_let == "A", 1, cap_join$Session )
cap_join$Session <-ifelse( cap_join$Session_let == "B", 2, cap_join$Session )
cap_join$Session <-ifelse( cap_join$Session_let == "C", 3, cap_join$Session )
cap_join$Session <-ifelse( cap_join$Session_let == "D", 4, cap_join$Session )

##remove Session_let and Site.Sess
cap_join <- cap_join[,c( "Individual", "Name", "Session", "East", "North")]




##################################################################
##### MERGE DATA
##################################################################

##merge site and capture data
FINAL<-merge(site_join, cap_join, by.x=c("SITE_NAME", "Session", "Easting.DC","Northing.DC"),  by.y=c("Name", "Session", "East", "North"), all=TRUE)


##Add in SEX
FINAL<-merge(FINAL, sex, by="Individual", all=TRUE)

###checks
length(unique(FINAL$SITE_NAME))==length(unique(site_join$SITE_NAME))
length(unique(FINAL[!is.na(FINAL$Individual),]$Individual))==length(unique(cap$Individual))
nrow(FINAL[!is.na(FINAL$Individual),])==nrow(cap)

##Add Trap Type
FINAL$Trap.Type<-"BS"

##Add Project
FINAL$Project<-"S_Selk_2005"

##Add Zone
FINAL$Zone<-11

##Add Year
FINAL$Year<-"2005"

##Create unique Site ID for each site
FINAL$SITE_NAME<-paste(FINAL$SITE_NAME,FINAL$Easting.DC,FINAL$Northing.DC,sep="_")


##Reorder
FINAL <- FINAL[, c("Project","Year", "SITE_NAME", "Trap.Type",  "Session",     "Easting.DC",  "Northing.DC", "Zone", "Individual"  , "Sex"   )]
FINAL<-FINAL[with(FINAL, order(SITE_NAME, Session)), ]

##EXPORT
write.csv(FINAL, "/Users/clayton.lamb/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/ProjectFiles/Single_Year/S_Selkirks_2005.csv", row.names = FALSE, na="")



####
unique(FINAL$SITE_NAME[!FINAL$SITE_NAME %in% bb$Site.Name])
as.character(FINAL$SITE_NAME[272])
