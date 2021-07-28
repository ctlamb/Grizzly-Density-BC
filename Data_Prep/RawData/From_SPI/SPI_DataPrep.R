#################################
## PREP SPI DATA FOR USE IN DENSITY META-ANALYSIS
#################################

##CLEAR WORKSPACE
rm(list = ls())

##Load Spatial packages
library(sf)
library(sp)
library(plyr)
library(ggmap)
library(mapview)
library(lubridate)
library(here)
library(tidyverse)


##################################################################
####LOAD CAPTURE DATA
##################################################################
cap <- ldply(list.files(path=here::here("inventories"), pattern="*.csv", full.names=TRUE), read.csv)


##################################################################
##### PLOT MAP ALL OF ALL SITE DATA
##################################################################

##Cull any missing spatial data
 cap<-subset(cap,!is.na(Easting.DC))

##one errant location for lower bowron (wrong UTM)
cap[cap$Study.Area.Name%in%"Lower Bowron River" & cap$Design.Component.Label %in% "s34_2B_U", "UTM.Zone.DC"] <- 10

##rename "Kingcome 1997"  to "Kingcome"  
cap <- cap%>%mutate(Study.Area.Name=as.character(Study.Area.Name))%>%
  mutate(Study.Area.Name=case_when(Study.Area.Name%in% "Kingcome 1997"~  "Kingcome", TRUE~Study.Area.Name))


##Loop to get each zone into correct UTM's and then project to lat long
Zones<-unique(cap$UTM.Zone.DC)
for(i in 1:length(Zones)){
  a<-subset(cap,UTM.Zone.DC==Zones[i])
  coordinates(a) <- c("Easting.DC","Northing.DC")
  proj4string(a)<-paste("+proj=utm +zone=",Zones[i]," +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",sep="")
  # spTransform function provide transformation between projections.
  a<- spTransform(a, CRS("+proj=utm +zone=11, +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  assign(paste("Cap",Zones[[i]],sep="_"), a)
}

##Join all CI Zones into a single layer
cap<-rbind(Cap_8,Cap_9,Cap_11,Cap_10)

##retain utm 11
cap$UTM.Zone.DC<-11
cap$Easting <- coordinates(cap)[,1]
cap$Northing <- coordinates(cap)[,2]

##transform to lat long
cap<-spTransform(cap,CRS("+init=epsg:4326"))
#mapview(cap["Study.Area.Name"])
cap<-as.data.frame(cap)

###Plot to make sure all is good
register_google("AIzaSyCOwGx2D77XOqRgGhKmcb5F4Kt_S61tCLI")
basemap <- get_map(location = c(lon = min(cap$Easting.DC)+(max(cap$Easting.DC)-min(cap$Easting.DC))/2, lat = min(cap$Northing.DC)+(max(cap$Northing.DC)-min(cap$Northing.DC))/2), zoom = 5)
map <- ggmap(basemap)+
  geom_point(aes(x = Easting.DC,
                 y = Northing.DC,
                 fill=Study.Area.Name),
             colour = "black",
             pch = 21,
             data=cap,size=5)+
  theme(legend.position = "none")
map

##################################################################
##### SET UP LOOP
##################################################################

Proj<-as.character(unique(cap$Study.Area.Name))

for(i in 1:length(Proj)){
  
  
  ##################################################################
  ##### PREP CAPTURE DATA
  ##################################################################
  
  ##Subset out project of interest
a <- cap%>%filter(Study.Area.Name%in%Proj[i])

##add year
a$Year <- a$Date%>%ymd()%>%year()

##drop w/missing year
a<- a%>%
  drop_na(Year)

  
  ###Plot to make sure all is good
basemap <- get_map(location = c(lon = median(a$Easting.DC), lat = median(a$Northing.DC)), zoom = 8)
map <- ggmap(basemap)+
  geom_point(aes(x = Easting.DC,
                 y = Northing.DC,
                 fill=Study.Area.Name),
             colour = "black",
             pch = 21,
             data=a,size=3)

print(map)


  ###Extract Sessiom
  a$Session <- str_sub(a$StudyAreaVisitNote, start=-1, end=-1)
  
  ##define name
  a$Sampling.Project=paste(Proj[i], a$Year, sep="_")
  
  ##Remove Opportunistic
  a <- subset(a,!Sample.Is.Incidental%in%"Y")

  ##retain useful columns
  a <- a%>%select("Sampling.Project", "Year", "Site.Design.Component.Label", "Trap.Type", "Session","Easting","Northing","UTM.Zone.DC", "Animal.ID", "Sex","Species")
  

  ##Remove factors
  a$Site.Design.Component.Label <- factor(a$Site.Design.Component.Label)
  a$Session <- factor(a$Session)
  a$Animal.ID<- factor(a$Animal.ID)
  
  ##black bear
bb <- a%>%mutate(bbear=case_when(Species%in%"M-URAM"~1, TRUE~0))%>%
  group_by(Sampling.Project, Year, Site.Design.Component.Label, Trap.Type, Session, Easting, Northing, UTM.Zone.DC)%>%
        summarise(bb.count=sum(bbear))
write_csv(bb,paste("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/CleanData/ProjectFiles/blackbear/",Proj[i],".csv", sep=""))
  
  
  #####Break into detect and nodetect data
  detect<-a%>%filter(!Animal.ID%in%"")%>%filter(Species%in% "M-URAR")
  NONdetect<-a%>%filter(Animal.ID%in%""| is.na(Animal.ID))
  
  ##Pull out unique records, based on certain columns while retaining all columns
  detect<-detect%>%distinct(Year, Site.Design.Component.Label, Session, Trap.Type, Animal.ID,
                            .keep_all = TRUE)
  
  ##remove NON detect site-sessions where a bear was detected
  detect$S_Sess <- paste(detect$Year, detect$Site.Design.Component.Label, detect$Session, sep="_")
  NONdetect$S_Sess <- paste(NONdetect$Year,NONdetect$Site.Design.Component.Label, NONdetect$Session, sep="_")
  NONdetect <- NONdetect[!NONdetect$S_Sess %in% detect$S_Sess,]

  ##Pull out unique records, based on certain columns while retaining all columns
  NONdetect<-NONdetect%>%distinct(Year, Site.Design.Component.Label, Session, Trap.Type, Animal.ID,
                                  .keep_all = TRUE)
  
  NONdetect$Animal.ID <- ""
  
  ##Make FINAL data
  FINAL<-rbind(detect,  NONdetect)
  

##Rename
FINAL<-FINAL%>%dplyr::rename("Project"="Sampling.Project",
                             "SITE_NAME"="Site.Design.Component.Label", 
                             "Zone"="UTM.Zone.DC",
                             "Individual"="Animal.ID",
                             "Easting.DC"="Easting", 
                             "Northing.DC"="Northing")


##Pull out unique records, based on certain columns while retaining all columns
FINAL<-FINAL%>%distinct(Year, SITE_NAME, Session, Trap.Type, Individual,
                     .keep_all = TRUE)



FINAL <- FINAL[, c("Project", "Year", "SITE_NAME", "Trap.Type",  "Session",     "Easting.DC",  "Northing.DC", "Zone",  "Individual"  , "Sex"   )]


##reclass trap type
FINAL<- FINAL%>%
  mutate(Trap.Type=case_when(Trap.Type%in%"HR-BA-BS"~"BS", !Trap.Type%in%"HR-BA-BS"~"RT"))


if( length(unique(FINAL$Year))==1){
  ##Reorder
FINAL<-FINAL[with(FINAL, order(SITE_NAME, Session)), ]

##EXPORT
write.csv(FINAL, paste("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/CleanData/ProjectFiles/Single_Year/",Proj[i], ".csv", sep=""), row.names = FALSE, na="")
}



if( length(unique(FINAL$Year))>1){
  ##Reorder
  FINAL<-FINAL[with(FINAL, order(Year, SITE_NAME, Session)), ]
  
  ##EXPORT
  write.csv(FINAL, paste("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/CleanData/ProjectFiles/Multi_Year/",Proj[i],"_", min(FINAL$Year), "_", max(FINAL$Year), ".csv", sep=""), row.names = FALSE, na="")
}
print(Proj[i])
}








