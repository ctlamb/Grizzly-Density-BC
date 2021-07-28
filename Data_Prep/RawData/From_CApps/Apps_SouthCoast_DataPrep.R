#################################
## PREP APPS CAPTURE DATA FOR USE IN DENSITY META-ANALYSIS
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

##################################################################
####LOAD CAPTURE DATA
##################################################################

#set working directory
setwd("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/RawData/From_CApps")


cap <- read.csv("Combined Detection Database_Apps_CLambClean.csv", header=TRUE, sep=",")


##################################################################
##### PLOT MAP ALL OF ALL SITE DATA
##################################################################

##Cull any missing spatial data
 cap<-subset(cap,!is.na(Easting))

####Make Data Spatial
coordinates(cap) <- c("Easting","Northing")
proj4string(cap)<-"+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"




##make sure UTMs aren't lost
cap$Easting.DC<-coordinates(cap)[,1]
cap$Northing.DC<-coordinates(cap)[,2]

##transform to lat long
cap<-spTransform(cap,CRS("+init=epsg:4326"))
cap<-as.data.frame(cap)

###Plot to make sure all is good

basemap <- get_map(location = c(lon = median(cap$Easting), lat = median(cap$Northing)), zoom = 7)
map <- ggmap(basemap)+
  geom_point(aes(x = Easting,
                 y = Northing,
                 fill=Sampling.Project),
             colour = "black",
             pch = 21,
             data=cap,size=5)

map


##################################################################
##### SELECT ONLY PROJECT OF INTEREST FROM SITE DATA AND LOOP

##FROM APPS:

#####The “monitoring” projects have some issues and are not appropriate for population estimation.

##USE
#Region-2004
#Region-2005
#Region-2006
#Region-2007
#Region-2007
#Toba-2008 
#Southgate-2010 


##CULL
#Monitoring-2005
#Monitoring-2006
#Monitoring-2007


##################################################################

##Rename Studies for Consistency with Other Projects
cap$Sampling.Project<-gsub("-", '_', cap$Sampling.Project, fixed = T)

##Extract Year
cap$Year <- sapply(strsplit(cap$Sampling.Project, "_"),`[`, 2)


##Drop studies from site data that are not good for abundance
cap<-subset(cap, !Sampling.Project %in% c("Monitoring_2005", "Monitoring_2006","Monitoring_2007"))



###Plot to make sure all is good


basemap <- get_map(location = c(lon = median(cap$Easting), lat = median(cap$Northing)), zoom = 7)
map <- ggmap(basemap)+
  geom_point(aes(x = Easting,
                 y = Northing,
                 fill=Sampling.Project),
             colour = "black",
             pch = 21,
             data=cap,size=3)

map




##################################################################
##### SET UP LOOP
##################################################################

Proj<-as.character(unique(cap$Sampling.Project))

for(i in 1:length(Proj)){
  
  
  ##################################################################
  ##### PREP CAPTURE DATA
  ##################################################################
  
  ##Subset out project of interest
  a <- cap[cap$Sampling.Project==Proj[i],]
  
  
  
  
  ###Plot to make sure all is good
  
  basemap <- get_map(location = c(lon = median(cap$Easting), lat = median(cap$Northing)), zoom = 7)
map <- ggmap(basemap)+
  geom_point(aes(x = Easting,
                 y = Northing,
                 fill=Sampling.Project),
             colour = "black",
             pch = 21,
             data=a,size=3)

print(map)

ggsave(paste("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/Maps/",Proj[i],".pdf",sep=""))


  ##retain useful columns
 
  a <- cap[cap$Sampling.Project==Proj[i], c("Year", "Trap.ID.per.site","Sessions","Easting.DC","Northing.DC", "Sess.1", "Sess.2", "Sess.3", "Sess.4","Sess.5")]
  
  ##Remove factors
  a$Trap.ID.per.site <- factor(a$Trap.ID.per.site)
  a$Sessions<- factor(a$Sessions)
  a$Sess.1<- factor(a$Sess.1)
  a$Sess.2<- factor(a$Sess.2)
  a$Sess.3<- factor(a$Sess.3)
  a$Sess.4<- factor(a$Sess.4)
  a$Sess.5<- factor(a$Sess.5)
  

  ##Melt dataframe so that each Sex is it's own row.  Retain Study, Site and Cell_Site info for each row
  b<-melt(a, na.rm=FALSE, id.vars=c("Year","Trap.ID.per.site","Sessions","Easting.DC","Northing.DC"), measured=c("Sess.1", "Sess.2", "Sess.3", "Sess.4","Sess.5"))
  
  ##Extract Session #
  b$Session<-substr(b$variable,start=nchar(as.character(b$variable)), stop=nchar(as.character(b$variable)))
  
  ##Cull Sessions not actually samples
  b<-b[b$Session<=as.numeric(as.character(b$Sessions))[1],]
  
  ###Split up value column as >1 bear captured/row
  c<-strsplit(b$value, ",")
  
  
  bearcol<-c("Bear1","Bear2","Bear3","Bear4","Bear5","Bear6")
  
  ##counter of max columns
  max<-1
  
  for(k in 1:length(c)){
    d<-c[[k]]
   
  for(j in 1:length(d)){
    if(length(d)>max){max<-length(d)}
      b[ k,bearcol[j] ]<-d[j]
    }
    
    
  }
  
  #####Break into detect and nodetect data
  detect<-b[!is.na(b[,c("Bear1")]),]
  NONdetect<-b[is.na(b[,("Bear1"),]),]
  
  ###checks
  if( any(!is.na(NONdetect[,c((ncol(NONdetect)-max+1):ncol(NONdetect))])) )
     warning(paste(Proj[i], "Splitting Caps and Non-Caps Failed at Non-Cap",sep=" "))
  
  ##NEED to make bear # flexible, and remove duplicated NA when bear is caught
  
  ##Melt AGAIN
  colnames(detect)
  detect
  
  drops <- c("Sessions", "variable", "value")
  detect<-detect[ , !(names(detect) %in% drops)]
  detectmelt<-melt(detect, na.rm=TRUE, id.vars=c("Year","Trap.ID.per.site","Session","Easting.DC","Northing.DC"), measured=c((ncol(detect)-max+1):ncol(detect)))
  
  ##Rbind detect and non-detect
  NONdetect<-NONdetect[,c("Year","Trap.ID.per.site","Session","Easting.DC","Northing.DC","variable","value")]
  
  e<-rbind(detectmelt,NONdetect)
  
  ##slim down data
  e<-e[,c("Year","Trap.ID.per.site", "Session", "Easting.DC", "Northing.DC", "value")]
  
  ##throw error is any duplicates
  if (any(duplicated(e)))
    warning(paste(Proj[i], "Duplicates after rbind",sep=" "))
  
  ###Pull out Sex Data
  e$Sex<-substr(e$value,start=1,stop=1)
  e[ !e$sex %in% c("M","F"), "Sex"]<-NA

  ##Make FINAL data
  FINAL<-e
  
  
##Add Trap Type
FINAL$Trap.Type<-"BS"

##Add Project
FINAL$Project<-as.character(Proj[i])

##Add Zone
FINAL$Zone<-10

##Rename
FINAL<-rename(FINAL,c("value"="Individual", "Trap.ID.per.site"="SITE_NAME"))


##Reorder
FINAL <- FINAL[, c("Project","Year", "SITE_NAME", "Trap.Type",  "Session",     "Easting.DC",  "Northing.DC", "Zone", "Individual"  , "Sex"   )]
FINAL<-FINAL[with(FINAL, order(SITE_NAME, Session)), ]

##EXPORT
write.csv(FINAL, paste("/Users/user/Dropbox/Work/gm_clayton/Analyses/Provincial_SECR/Data/CleanData/ProjectFiles/Single_Year/",Proj[i], ".csv", sep=""), row.names = FALSE, na="")

}

