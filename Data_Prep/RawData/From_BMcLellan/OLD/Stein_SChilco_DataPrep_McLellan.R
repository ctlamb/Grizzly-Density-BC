#################################
## PREP MCLELLAN COAST DATA FOR USE IN DENSITY META-ANALYSIS
#################################

##CLEAR WORKSPACE
rm(list = ls())

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
library(reshape2)
require(plotKML)
require(RFigisGeo)

##################################################################
####LOAD CAPTURE DATA
##################################################################

#set working directory
setwd("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_BMcLellan")


cap <- read.csv("McLellan South Coast 2010 to 2014.csv", header=TRUE, sep=",")




##################################################################
##### ADD SITE DATA TO ALL RECORDS, AS BRUCE HAND BALMED NON-CAPTURE RECORDS
##################################################################

###Remove capture records for now, retain RT's and HT's
cap <- subset( cap,!HTorRT %in% c("C","C (W_HELP)", "KILL", "KILL "))

##fix a few site names that don't start with HT or RT
cap$Site.Name <-  as.character(cap$Site.Name)
cap[cap$Site.Name=="Tenquil Trap", "Site.Name"] <- "RT Tenquil Trap"
cap[cap$Site.Name=="W Texas trap site", "Site.Name"] <- "RT W Texas trap site"
cap[cap$Site.Name=="Whitebark Pine Tree in Moly", "Site.Name"] <- "RT Whitebark Pine Tree in Moly"


##fix a few site names that have year in them
cap[cap$Site.Name=="HT Hope Slide2010", "Site.Name"] <- "HT Hope Slide"
cap[cap$Site.Name=="HT Hope Slide2012", "Site.Name"] <- "HT Hope Slide"


cap[cap$Site.Name=="RT McGill N 2010_2", "Site.Name"] <- "RT McGill N 1"

cap[cap$Site.Name=="RT  U Hope ", "Site.Name"] <- "RT Hope1"

cap[cap$Site.Name=="RT E Gott1", "Site.Name"] <- "RT E Gott"

summary(as.factor(cap$Site.Name))

##
##extract site type
##
cap$Trap.Type <- substr(cap$Site.Name, start=0, stop=2)

##make all site types the same
cap$Trap.Type <- ifelse(cap$Trap.Type %in% c("HT","Ht"),"BS",as.character(cap$Trap.Type))
cap$Trap.Type <- ifelse(cap$Trap.Type %in% c("RT","Rt"),"RT",as.character(cap$Trap.Type))

##check
summary(as.factor(cap$Trap.Type ))

###Remove empty records
cap <- subset( cap,Trap.Type != "" )

##check
summary(as.factor(cap$Trap.Type ))

cap$Trap.Type <- as.factor(cap$Trap.Type )


##################################################################
##### CLEAN UP IDS
##################################################################
cap$NEWIDS <- as.character(cap$NEWIDS)
cap$ID <- as.character(cap$ID)

cap[cap$NEWIDS=="NONE","NEWIDS"] <- ""

cap[cap$NEWIDS!="" & cap$ID=="", "ID"] <- cap[cap$NEWIDS!="" & cap$ID=="", "NEWIDS"]

cap[cap$NEWIDS!="" & cap$ID==" ", "ID"] <- cap[cap$NEWIDS!="" & cap$ID==" ", "NEWIDS"]


cap[cap$ID=="M320", "ID"]
cap[cap$ID==" M320", "ID"]


################################################
##Spatially assign GBPU
################################################

##Cull any missing spatial data
cap <- subset(cap,!is.na(Easting))


####Make Data Spatial
coordinates(cap) <- c("Easting","Northing")
proj4string(cap) <- "+proj=utm +zone=10+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"


##make sure UTMs aren't lost
cap$Easting.DC<-coordinates(cap)[,1]
cap$Northing.DC<-coordinates(cap)[,2]

##Load GBPU
GBPU <- readOGR(dsn="/Users/user/Documents/University/Geographic_Data/Bears/GBPU_Boundaries_Bears/GBPU_BC", "GBPU_BC_Only")
GBPU <- spTransform(GBPU, CRS(proj4string(cap))) 

##plot
plot(GBPU)
plot(cap,add=TRUE,col="red")

##assign GBPU
cap$GBPU <- over(cap, GBPU)$GBPU_NAME


##Shorten
cap$GBPU <- ifelse(cap$GBPU=="South Chilcotin Ranges","SC",as.character(cap$GBPU))
cap$GBPU <- ifelse(cap$GBPU=="Stein-Nahatlatch","SN",as.character(cap$GBPU))


##################################################################
##### PREP STUDY BOUNDS
##################################################################

##Loading .KMLs

##SC
SC_SA<-readOGR(dsn="/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_BMcLellan/Stein_Spatial/SouthChilco_DNAsampling.kml", "SouthChilco_DNAsampling.kml")

##SN
SN_SA<-readOGR(dsn="/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_BMcLellan/Stein_Spatial/SteinNahatl_DNAsampling.kml", "SteinNahatl_DNAsampling.kml")

##merge poly to join SC & SN together
merge<-readOGR(dsn="/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_BMcLellan/Stein_Spatial/Coast_MergePoly.kml", "Coast_MergePoly.kml")



##Merge all together
SA <- gUnion(gUnion(SC_SA,merge),SN_SA)


##reproject
SA <- spTransform(SA, CRS(proj4string(cap))) 

###CLIP Polygon by Polygon to get SA with GBPU polys
GBPU_SC_SN <- subset(GBPU, GBPU_NAME %in% c("Stein-Nahatlatch", "South Chilcotin Ranges") & GBPU_VERS == 2012)
SA_GBPU <- intersect(SA, GBPU_SC_SN)


##Write to .shp
writeOGR(SA_GBPU, "/Users/user/Documents/University/U_A/Analyses/BM_Coast_Pradel/Coast_SECR/BMCoast_SA", "coastSA_w_GBPU", driver="ESRI Shapefile",overwrite=TRUE)


##fortify for ggplot2
SA_GBPU_fort<- fortify(spTransform(SA_GBPU,CRS("+init=epsg:4326")))


##################################################################
##### PLOT MAP ALL OF ALL SITE DATA
##################################################################


##transform to lat long
cap<-spTransform(cap,CRS("+init=epsg:4326"))
cap<-as.data.frame(cap)



###Plot to make sure all is good

basemap <- get_map(location = c(lon = median(cap$Easting), lat = median(cap$Northing)), zoom = 9)
map <- ggmap(basemap)+
  geom_point(aes(x = long,
                 y = lat),
             fill=NA,
             data=SA_GBPU_fort)+
  geom_point(aes(x = Easting,
                 y = Northing,
                 fill=GBPU),
             colour = "black",
             pch = 21,
             data=cap,size=4)

map


map <- ggmap(basemap)+
  geom_point(aes(x = Easting,
                 y = Northing),
             colour = "black",
             fill="red",
             pch = 21,
             data=cap,size=1)+
         facet_grid(. ~ YR)
  

#map

#ggsave("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_BMcLellan/Plots_Interm_DELETE/SC_SN_ByYear.pdf", width=12, height=4,units="in")

###Remove years of data that Bruce says are just Rochertta's "Monitoring"
##Good years are:
# So, SN, 2005, 2010, 2011, 2012, 2013, 2014 
# 
# and SC 2006, 2010, 2011, 2012, 2013, 2014


  
  ##################################################################
  ##### PREP CAPTURE DATA
  ##################################################################
 
##Project name
cap$Project <- "Lilloet_2005_2015"

##UTM zone
cap$UTM.Zone.DC <- 10

##UTM zone
cap$UTM.Zone.DC <- 10


##retain useful columns
a <- cap[, c("Project", "Site.Name", "Trap.Type", "Session","Easting.DC","Northing.DC","UTM.Zone.DC", "ID", "HIT","YR", "Sex")]






  ##Remove factors
  a$Site.Name <- factor(a$Site.Name)
  a$Session <- factor(a$Session)
  a$ID<- factor(a$ID)

  #####Break into detect and nodetect data
  detect<-a[a$ID!="",]
  NONdetect<-a[a$ID=="",]
  
  ##Pull out unique records, based on certain columns while retaining all columns
  detect<-detect[!duplicated(detect[,c("YR","Site.Name", "Session", "Trap.Type", "ID")]),]
  
  ##remove NON detect site-sessions where a bear was detected
  detect$S_Sess <- paste(detect$Site.Name, detect$Session, detect$YR, sep="_")
  NONdetect$S_Sess <- paste(NONdetect$Site.Name, NONdetect$Session, NONdetect$YR, sep="_")
  NONdetect <- NONdetect[!NONdetect$S_Sess %in% detect$S_Sess,]

  ##Pull out unique records, based on certain columns while retaining all columns
  NONdetect<-NONdetect[!duplicated(NONdetect[,c("YR", "Site.Name", "Session", "Trap.Type", "ID")]),]
  
  ##Make FINAL data
  FINAL<-rbind(detect,  NONdetect)
  


##want
##Project	SITE_NAME	Trap.Type	Session	Easting.DC	Northing.DC	Zone	Individual	Sex

##Rename
FINAL<-rename(FINAL,c("Site.Name"="SITE_NAME", 
                      "UTM.Zone.DC"="Zone",
                      "YR"="Year",
                      "ID" = "Individual"))
                      

##Reorder
FINAL <- FINAL[, c("Project", "Year", "SITE_NAME", "Trap.Type",  "Session",     "Easting.DC",  "Northing.DC", "Zone", "Individual"  , "Sex"   )]
FINAL<-FINAL[with(FINAL, order(SITE_NAME, Session)), ]


##Load SN and SC 2005-2006
SN05 <- read.csv("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/ProjectFiles/Single_Year/Region_2005.csv", header=TRUE, sep=",")
SC05 <- read.csv("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/ProjectFiles/Single_Year/Region_2006.csv", header=TRUE, sep=",")

##Bind onto FINAL
FINAL <- rbind(FINAL, rbind(SN05,SC05))

##Project name
FINAL$Project <- "Lilloet_2005_2015"

####Make Data Spatial
coordinates(FINAL) <- c("Easting.DC","Northing.DC")
proj4string(FINAL) <- "+proj=utm +zone=10+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

##clip to McLellan SA only
FINAL$InSA <- over(FINAL, SA)
FINAL_cull<-subset(FINAL,!is.na(InSA))


##remove inSA and turn back to data.frame
FINAL_cull <- as.data.frame(FINAL_cull)
FINAL_cull <- subset(FINAL_cull, select=-InSA)

##Add in missing Sex  
FINAL_cull[FINAL_cull$Sex=="" & FINAL_cull$Individual!="", "Sex"] <- substr(FINAL_cull[FINAL_cull$Sex=="" & FINAL_cull$Individual!="", "Individual"], start = 0, stop =1)
FINAL_cull[FINAL_cull$Sex=="" & FINAL_cull$Individual!="", "Sex"] <- substr(FINAL_cull[FINAL_cull$Sex=="" & FINAL_cull$Individual!="", "Individual"], start = 0, stop =1)




##EXPORT
write.csv(FINAL_cull, "/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/ProjectFiles/Multi_Year/Lilloet_2005_2014.csv", row.names = FALSE, na="")










