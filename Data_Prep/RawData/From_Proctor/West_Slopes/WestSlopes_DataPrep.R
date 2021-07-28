#################################
## PREP PROCTOR CAPTURE DATA FOR USE IN DENSITY META-ANALYSIS
## WEST-SLOPES (GOLDEN)
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
require(magrittr)
require(dplyr)

##################################################################
####LOAD CAPTURE DATA
##################################################################

#set working directory
setwd("/Users/clayton.lamb/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_Proctor/West_Slopes")


cap <- read.csv("WS all genos LATEST BEST_CLEAN_LAMB.csv", header=TRUE, sep=",")
cap_sex <-  read.csv("WS filter from MP Master genotype file PhD_CLEAN_LAMB.csv", header=TRUE, sep=",")
site <- read.csv("WS site-session data 96-98_apps.csv", header=TRUE, sep=",") 


##################################################################
##### PLOT MAP ALL OF ALL SITE DATA
##################################################################

##Cull any missing spatial data
site<-subset(site,!is.na(easting))

####Make Data Spatial
coordinates(site) <- c("easting","northing")
proj4string(site)<-"+proj=utm +zone=11 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"



##Cull any missing data

site$Easting.DC<-coordinates(site)[,1]
site$Northing.DC<-coordinates(site)[,2]

##transform to lat long
site<-spTransform(site,CRS("+init=epsg:4326"))
site<-as.data.frame(site)

###Plot to make sure all is good
site$Year <- as.factor(site$Year)

basemap <- get_map(location = c(lon = median(site$easting), lat = median(site$northing)), zoom = 8)
map <- ggmap(basemap)+
  geom_point(aes(x = easting,
                 y = northing,
                 fill=Year),
             colour = "black",
             pch = 21,
             data=site,size=5)

map



##################################################################
##### CLEAN UP CAP DATA
##################################################################


##Remove cap records with no individual
cap <- subset(cap, Perm.ID!="")

##Remove any individuals w/ question marks:
matches <- 
cap <- cap[!grepl("\\?", cap$Perm.ID),]

##Pull out unique records, based on certain columns while retaining all columns
Ucap<-cap[!duplicated(cap[,c("Perm.ID","Year", "Site" ,"Session" )]),]

##retain useful columns
Ucap <- Ucap[,c("Perm.ID","Year", "Site" ,"Session")]

##how many bears?
length(unique(cap$Perm.ID))

##clean sex
Ucap_sex <- cap[,c("Perm.ID","Sex")] %>%
            rbind(cap_sex[,c("Perm.ID","Sex")]) %>%
            subset(Sex %in% c("M","F")) %>%
            subset(!Perm.ID %in% c("315", "463")) %>% ##Bears =315 & 463 have 2 sexes
            unique()

##Join Sex
cap_sexjoin <- merge(Ucap, Ucap_sex, by="Perm.ID", all=TRUE)

##So, who doesn't have sex ID'd?
Ucap$Perm.ID [!Ucap$Perm.ID %in% Ucap_sex$Perm.ID]


##Cull non-cap records and fill missing sex w/U
cap_sexjoin <- subset(cap_sexjoin,!is.na(Site))
cap_sexjoin[is.na(cap_sexjoin$Sex),"Sex"] <- "U"


  
  ##################################################################
  ##### PREP CAPTURE DATA
  ##################################################################
  
##Join Capture and Site Data
site$Site <- as.factor(site$Site)
site$Session <- as.factor(site$Session)
site$Year <- as.factor(paste("19",site$Year,sep=""))

summary(site$Year)

cap_sexjoin$Site <- as.factor(cap_sexjoin$Site)
cap_sexjoin$Session <- as.factor(cap_sexjoin$Session)
cap_sexjoin$Year <- as.factor(cap_sexjoin$Year)

colnames(site)
cap_join <- merge(site[,c("Site","Session", "Year","Easting.DC",  "Northing.DC")], cap_sexjoin, by=c("Site", "Session","Year"), all=TRUE)

length(unique(cap_join$Perm.ID))-1
length(unique(cap$Perm.ID))


length(unique(cap_join$Easting.DC))
length(unique(site$Easting.DC))
  

###Change sessions from letters to #'s, A's are equlivelent to sess 1, B=2, C=3, D=4
  cap_join$Session<-as.character(cap_join$Session)
  cap_join$Session <-ifelse( cap_join$Session == "A", as.character(1), cap_join$Session )
  cap_join$Session <-ifelse( cap_join$Session == "B", as.character(2), cap_join$Session )
  cap_join$Session <-ifelse( cap_join$Session == "C", as.character(3), cap_join$Session )
  cap_join$Session <-ifelse( cap_join$Session == "D", as.character(4), cap_join$Session )
  
  
##change df name
  FINAL <-  cap_join

##Add Trap Type
FINAL$Trap.Type<-"BS"

##Add Project
FINAL$Project<-"West_Slopes_1996_1998"

##Add Zone
FINAL$Zone<-11

##Rename
FINAL<-plyr::rename(FINAL,c("Perm.ID"="Individual"))

##Create unique Site ID for each site
FINAL$SITE_NAME<-paste(FINAL$Site,FINAL$Session,FINAL$Year,sep="_")

##Reorder
FINAL <- FINAL[, c("Project","Year", "SITE_NAME", "Trap.Type",  "Session",     "Easting.DC",  "Northing.DC", "Zone", "Individual"  , "Sex"   )]
FINAL<-FINAL[with(FINAL, order(SITE_NAME, Session)), ]

##EXPORT
write.csv(FINAL, "/Users/clayton.lamb/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/ProjectFiles/Multi_Year/West_Slopes_1996_1998.csv", row.names = FALSE, na="")





##################################################################
##### SUMMARY STATS
##################################################################
subset(FINAL, !is.na(Individual))%>%
group_by(Year) %>%
  summarise(Unique_Elements = n_distinct(Individual))