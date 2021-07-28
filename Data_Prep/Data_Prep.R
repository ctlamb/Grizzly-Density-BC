#################################
## LOAD ALL DNA DATA
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
library(data.table)
library(ggsn)
library(ggplot2)
library(here)
library(adehabitatHR)
library(tidyverse)

##################################################################
####LOAD CAPTURE DATA
##################################################################

temp1 = list.files(path=here::here("CleanData", "ProjectFiles","Single_Year"),pattern = "*.csv", full.names=TRUE)
temp2 = list.files(path=here::here("CleanData", "ProjectFiles","Multi_Year"),pattern = "*.csv", full.names=TRUE)

temp <- append(temp1, temp2)


df <- do.call(rbind, lapply(temp, fread))


##change individual NA's to ""
df<- df%>%mutate(Individual=case_when(is.na(Individual)~"", TRUE~Individual))
##################################################################
##### PLOT MAP ALL OF ALL SITE DATA
##################################################################

##Cull any missing spatial data
df<-subset(df,!is.na(Easting.DC))


Zones<-unique(df$Zone)


##Loop to get each zone into correct UTM's and then project to lat long
for(i in 1:length(Zones)){
  a<-subset(df,Zone==Zones[i])
  coordinates(a) <- c("Easting.DC","Northing.DC")
  proj4string(a)<-paste("+proj=utm +zone=",Zones[i]," +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",sep="")
  # spTransform function provide transformation between projections.
  a<- spTransform(a, CRS("+proj=longlat +ellps=WGS84"))
  assign(paste("Cap",Zones[[i]],sep="_"), a)
}

##Join all CI Zones into a single layer
All<-rbind(Cap_11,Cap_10)

##back to DF
All<-as.data.frame(All)

###Plot to make sure all is good
register_google("AIzaSyCOwGx2D77XOqRgGhKmcb5F4Kt_S61tCLI")
basemap <- get_map(location = c(lon = min(All$Easting.DC)+(max(All$Easting.DC)-min(All$Easting.DC))/2, lat = min(All$Northing.DC)+(max(All$Northing.DC)-min(All$Northing.DC))/2), zoom = 5)
map <- ggmap(basemap)+
  geom_point(aes(x = Easting.DC,
                 y = Northing.DC,
                 fill=Project),
             colour = "black",
             pch = 21,
             data=All,size=2)+
  theme(legend.position = "none")

map
ggsave(here::here("CleanData", "Maps", "all.png"), width=10, height=10, units="in")



##################################################################
### SUMMARY STATS
##################################################################

##N Projects
length(unique(df$Project))


##N individuals
length(unique(df$Individual))-1

##N captures
sum(df$Individual!="")


##N Sites
length(unique(df$SITE_NAME))





##################################################################
### REPROJECT INTO COMMON COORDINATE SYSTEM
##################################################################
##Join all CI Zones into a single layer
All2<-rbind(Cap_11,Cap_10)
plot(All2)


###TRANSFORM INTO BC ALBERS


##Code: EPSG::3005
##Name: NAD83 / BC Albers
##From: http://www.epsg-registry.org

##from BC site:
#BC Environment GIS Working Group has chosen a standard projection and datum for all spatial data stored in ARC/INFO.
#The projection is Albers Equal Area Conic, with parameters of :

#  Central meridian: -126.0 (126:00:00 West longitude)
#First standard parallel: 50.0 (50:00:00 North latitude)
#Second standard parallel: 58.5 (58:30:00 North latitude)
#Latitude of projection origin: 45.0 (45:00:00 North latitude)
#False northing: 0.0
#False easting: 1000000.0 (one million metres)
#The datum is NAD83, based on the GRS80 ellipsoid.

#Its distance representation is also very close to Universal Transverse Mercator. 
#We found the average difference in latitude to be 0.15 %, and the maximum difference to be 0.30 %
#within the province. The average difference in longitude is 0.17 %, and the maximum difference is 0.32 %.


BC_ALL<- spTransform(All2, CRS("+init=epsg:3005"))

plot(BC_ALL, col="red")




##################################################################
### DO SOME DATA MANAGEMENT
##################################################################

##standardize U and UC sex
BC_ALL$Sex <- ifelse(BC_ALL$Sex%in%c("UC", "U"),NA,as.character(BC_ALL$Sex ))

##Convert Individual, Sex, SITE_NAME, and Trap.Type to factors
BC_ALL$Individual <-as.character(BC_ALL$Individual)
BC_ALL$Sex <-as.factor(BC_ALL$Sex)
BC_ALL$Trap.Type <-as.factor(BC_ALL$Trap.Type)

##Convert SITE_NAME to universal name (X_Y) then to Factor
BC_ALL$SITE_NAME <- paste(abs(round(BC_ALL$Easting.DC,0)),round(BC_ALL$Northing.DC,0), BC_ALL$Trap.Type, sep="_")
BC_ALL$SITE_NAME <-as.factor(BC_ALL$SITE_NAME)

##rename sites, as SECR doesn't like current names
#BC_ALL <- transform(BC_ALL, SITE_NAME=as.factor(match(SITE_NAME, unique(SITE_NAME))))

##remove whitepace from project names
BC_ALL$Project <- str_replace_all(BC_ALL$Project, " ", "_")

##remove YEAR
BC_ALL$Project.group <- str_sub(BC_ALL$Project, start=0,end=-6)


##address naming issue, some bears called same name between projects, but not actually same bear
BC_ALL$Individual <- ifelse(!is.na(BC_ALL$Individual) & BC_ALL$Individual!="",  paste(gsub('[0-9]+', '',sapply(strsplit(BC_ALL$Project.group, "_"), function(x){
  paste(substring(x, 1, 1), collapse = "")})),  BC_ALL$Individual, sep="_"),
  as.character(BC_ALL$Individual))





##################################################################
### SUMMARY STATS
##################################################################

##N Projects
length(unique(BC_ALL$Project))
length(unique(BC_ALL$Project.group))

##N individuals
length(unique(BC_ALL$Individual))-1

##N captures
sum(BC_ALL$Individual!="")


##N Sites
length(unique(BC_ALL$SITE_NAME))











##################################################################
  ### LOOP THROUGH PROJECTS AND PLOT MOVEMENTS
##################################################################
Projs<-unique(BC_ALL$Project)



##Transform into GEO
BC_ALL_geo <- spTransform(BC_ALL, CRS("+init=epsg:4326"))

##LOOP
for(i in 1:length(Projs)){
  a<-as.data.frame(BC_ALL_geo[BC_ALL_geo$Project==Projs[i],])

  
  ##make bounding box to extract google maps data with
bbox <- make_bbox(Easting.DC, Northing.DC, a, f = 0.8)

##then one to clip with
bbox_clip <- make_bbox(Easting.DC, Northing.DC, a, f = 0.2)

##set up scalar for point size
xrange <- max(a$Easting.DC)-min(a$Easting.DC)+1
yrange <- max(a$Northing.DC)-min(a$Northing.DC)+1

scalar<-(xrange*yrange)/2.255847


##prep info about project, sessions, captures, etc

occ<-paste( length(unique(paste(a$Year,a$Session,sep="_"))) ,"Occassions,",sep=" ")
det<-paste(nrow(a[a$Individual!="",]),"Detections,",sep=" ")
ind<-paste(length(unique(a[a$Individual!="",]$Individual)),"Individuals",sep=" ")

line<-paste(occ,det,ind,sep=" ")

##PLOT

basemap <- get_map(location = bbox)
map <- ggmap(basemap)+
  geom_point(aes(x = Easting.DC,
                 y = Northing.DC,
                 fill= "gold2"),
             colour = "black",
             pch = 21,
             data=a[a$Individual=="",],size=5-scalar+1)+
  geom_path(aes(x = Easting.DC,
                y = Northing.DC,
                group=Individual),
            colour = "black",
            data=a[a$Individual!="",],
            lwd=1.5)+
  geom_path(aes(x = Easting.DC,
                y = Northing.DC,
                group=Individual,
                colour = "firebrick2"),
            data=a[a$Individual!="",],
            lwd=1)+

  geom_point(aes(x = Easting.DC,
                 y = Northing.DC,
                 fill= "indianred1"),
             colour = "black",
             pch = 21,
             data=a[a$Individual!="",],size=6-scalar+1)+
  ggtitle(bquote(atop(.(Projs[i]), atop(italic(.(line)), "")))) +
  #ggtitle( paste0(paste0(Projs[i]),"\n", paste0(line)) )+
  theme(plot.title = element_text(size = 25, face = "bold", colour = "black", vjust = 1),
        legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
        legend.key = element_blank())+
  coord_map(xlim = c(bbox_clip[1], bbox_clip[3]),ylim = c(bbox_clip[2], bbox_clip[4]))+
  scale_fill_manual(name = 'Grizzly Detected',  guide = 'legend', values=c("gold2","firebrick2"),labels = c('No','Yes'))+
  scale_color_manual(name="", guide = 'legend', values=c("firebrick2"),labels = c("Recapture Path"))


ggsave(here::here("CleanData","Maps", "byproj", paste(Projs[i], "png", sep=".")), width=10, height=10, units="in")

}


##################################################################
### LOOP THROUGH PROJECTS AND PLOT STUDY BOUNDS
##################################################################

Projs<-unique(BC_ALL$Project)


shp <- list()
##LOOP
for(i in 1:length(Projs)){
  a<-BC_ALL[BC_ALL$Project==Projs[i],]
  
  ## Make MCP
  MCP<- mcp(a, percent=95)
  ##buffer a bit
  MCP<-gBuffer(MCP, width=5000)
  

  ##put into list
  shp[[i]] <-  MCP
    
}

##Name
names(shp) <- paste0(Projs)

##prep names
IDs <- names(shp)

##Combine to SpatialPolygons
Spol1 <- SpatialPolygons(lapply(1:length(shp), function(i) {
  Pol <- slot(shp[[i]], "polygons")[[1]]
  slot(Pol, "ID") <- IDs[i]
 Pol 
})
)

##To SpatialPolygonsDataFrame
SAs<-SpatialPolygonsDataFrame(Spol1,data.frame(Study=IDs),match.ID=FALSE)

##Set Projection
proj4string(SAs) <- proj4string(BC_ALL)
plot(SAs)
##export to .shp
writeOGR(SAs, here::here("Spatial_Layers_ProvSECR", "traps"), "DNA_StudyBounds", driver="ESRI Shapefile",overwrite=TRUE)



##################################################################
### PLOT ALL MOVEMENTS and BY REGION
##################################################################


Projs_Grouped<-list(unique(BC_ALL$Project.group),
                    c("Jumbo", "PWC", "S_Purcell", "HWY"),
                    c("Region","Southgate", "Toba", "Stein_Nahat"),
                    c("Southern_Rockies", "Flathead"))

Projs_Names <- c("British Columbia Grizzly DNA Projects",
                 "Purcells Grizzly DNA Projects",
                 "South Coast Grizzly DNA Projects",
                 "South Rockies Grizzly DNA Projects")


##Transform into GEO
BC_ALL_geo <- spTransform(BC_ALL, CRS("+init=epsg:4326"))

##LOOP
for(i in 1:length(Projs_Grouped)){
  
  
  a<-as.data.frame(BC_ALL_geo[BC_ALL_geo$Project.group %in% Projs_Grouped[[i]],])
  
  
  ##make bounding box to extract google maps data with
  bbox <- make_bbox(Easting.DC, Northing.DC, a, f = 0.8)
  
  ##then one to clip with
  bbox_clip <- make_bbox(Easting.DC, Northing.DC, a, f = 0.2)
  
  ##set up scalar for point size
  xrange <- max(a$Easting.DC)-min(a$Easting.DC)+1
  yrange <- max(a$Northing.DC)-min(a$Northing.DC)+1
  
  scalar<-sqrt(35.11358/(xrange*yrange))


  ##prep info about project, sessions, captures, etc
  
  occ<-paste( length(unique(paste(a$Project, a$Year,a$Session,sep="_"))) ,"Occassions,",sep=" ")
  det<-paste(nrow(a[a$Individual!="",]),"Detections,",sep=" ")
  ind<-paste(length(unique(a[a$Individual!="",]$Individual)),"Individuals",sep=" ")
  
  line<-paste(occ,det,ind,sep=" ")
  
  ##PLOT
  
  basemap <- get_map(location = bbox)
  map <- ggmap(basemap)+
    geom_point(aes(x = Easting.DC,
                   y = Northing.DC,
                   fill= "gold2"),
               colour = "black",
               pch = 21,
               data=a[a$Individual=="",],size=0.65*scalar)+
    geom_path(aes(x = Easting.DC,
                  y = Northing.DC,
                  group=Individual),
              colour = "black",
              data=a[a$Individual!="",],
              lwd=0.7*scalar)+
    geom_path(aes(x = Easting.DC,
                  y = Northing.DC,
                  group=Individual,
                  colour = "firebrick2"),
              data=a[a$Individual!="",],
              lwd=0.4*scalar)+
    
    geom_point(aes(x = Easting.DC,
                   y = Northing.DC,
                   fill= "indianred1"),
               colour = "black",
               pch = 21,
               data=a[a$Individual!="",],size=1.25*scalar)+
    ggtitle(bquote(atop(.(Projs_Names[i]), atop(italic(.(line)), "")))) +
    #ggtitle( paste0(paste0(Projs[i]),"\n", paste0(line)) )+
    theme(plot.title = element_text(size = 25, face = "bold", colour = "black", vjust = 1),
          legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
          legend.key = element_blank())+
    
    coord_map(xlim = c(bbox_clip[1], bbox_clip[3]),ylim = c(bbox_clip[2], bbox_clip[4]))+
    
    scale_fill_manual(name = 'Grizzly Detected',  guide = 'legend', values=c("gold2","firebrick2"),labels = c('No','Yes'))+
    
    scale_color_manual(name="", guide = 'legend', values=c("firebrick2"),labels = c("Recapture Path"))
  
  
  ggsave(here::here("CleanData","Maps", "grouped", paste(Projs_Names[i], "png", sep=".")), width=10, height=10, units="in")
  
  
}






##################################################################
### PLOT ALL MOVEMENTS and BY REGION FEMALES ONLY
##################################################################




Projs_Grouped<-list(unique(BC_ALL$Project.group),
                    c("Jumbo", "PWC", "S_Purcell", "HWY"),
                    c("Region","Southgate", "Toba", "Stein_Nahat"),
                    c("Southern_Rockies", "Flathead"))

Projs_Names <- c("British Columbia Grizzly DNA Projects",
                 "Purcells Grizzly DNA Projects",
                 "South Coast Grizzly DNA Projects",
                 "South Rockies Grizzly DNA Projects")



##LOOP
for(i in 1:length(Projs_Grouped)){
  
  
  a<-as.data.frame(BC_ALL_geo[BC_ALL_geo$Project.group %in% Projs_Grouped[[i]] & BC_ALL_geo$Sex %in% c("F", ""),])
  
  
  ##make bounding box to extract google maps data with
  bbox <- make_bbox(Easting.DC, Northing.DC, a, f = 0.8)
  
  ##then one to clip with
  bbox_clip <- make_bbox(Easting.DC, Northing.DC, a, f = 0.2)
  
  ##set up scalar for point size
  xrange <- max(a$Easting.DC)-min(a$Easting.DC)+1
  yrange <- max(a$Northing.DC)-min(a$Northing.DC)+1
  
  scalar<-sqrt(35.11358/(xrange*yrange))
  
  
  ##prep info about project, sessions, captures, etc
  
  occ<-paste( length(unique(paste(a$Project, a$Year,a$Session,sep="_"))) ,"Occassions,",sep=" ")
  det<-paste(nrow(a[a$Individual!="",]),"Detections,",sep=" ")
  ind<-paste(length(unique(a[a$Individual!="",]$Individual)),"Individuals",sep=" ")
  
  line<-paste(occ,det,ind,sep=" ")
  
  ##PLOT
  
  basemap <- get_map(location = bbox)
  map <- ggmap(basemap)+
    geom_point(aes(x = Easting.DC,
                   y = Northing.DC,
                   fill= "gold2"),
               colour = "black",
               pch = 21,
               data=a[a$Individual=="",],size=0.65*scalar)+
    geom_path(aes(x = Easting.DC,
                  y = Northing.DC,
                  group=Individual),
              colour = "black",
              data=a[a$Individual!="",],
              lwd=0.7*scalar)+
    geom_path(aes(x = Easting.DC,
                  y = Northing.DC,
                  group=Individual,
                  colour = "firebrick2"),
              data=a[a$Individual!="",],
              lwd=0.4*scalar)+
    
    geom_point(aes(x = Easting.DC,
                   y = Northing.DC,
                   fill= "indianred1"),
               colour = "black",
               pch = 21,
               data=a[a$Individual!="",],size=1.25*scalar)+
    ggtitle(bquote(atop(.(Projs_Names[i]), atop(italic(.(line)), "")))) +
    #ggtitle( paste0(paste0(Projs[i]),"\n", paste0(line)) )+
    theme(plot.title = element_text(size = 25, face = "bold", colour = "black", vjust = 1),
          legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
          legend.key = element_blank())+
    
    coord_map(xlim = c(bbox_clip[1], bbox_clip[3]),ylim = c(bbox_clip[2], bbox_clip[4]))+
    
    scale_fill_manual(name = 'Female Grizzly Detected',  guide = 'legend', values=c("gold2","firebrick2"),labels = c('No','Yes'))+
    
    scale_color_manual(name="", guide = 'legend', values=c("firebrick2"),labels = c("Recapture Path"))
  
  
  ggsave(here::here("CleanData","Maps", "grouped", paste("F",Projs_Names[i], "png", sep=".")), width=10, height=10, units="in")
  
  
}





########################################
### Now put into SECR format
########################################

Proj <- unique(BC_ALL$Project)
for(i in 1:length(Proj)){
  
df <- BC_ALL%>%as.data.frame()%>%filter(Project%in%Proj[i])


##Add in "Session Column" for SECR
df$`#Session`<-df$Project


##Retain columns needed for Capture Data,
CapData<-df[c("#Session","SITE_NAME","Session","Individual","Sex")]

##Remove all non-capture data
CapData<-subset(CapData, Individual!="")

##Rename Columns
CapData<-plyr::rename(CapData, c("Session"="Occassion", "SITE_NAME"="Detector", "Individual"="ID"))

##Reorder so the columns are Session, ID, Occassion and then Detector as per SECR format
CapData <- subset(CapData, select=c("#Session", "ID", "Occassion", "Detector","Sex"))

##remove sexes != M or F
#CapData<-subset(CapData,Sex %in% c("M","F"))

##Remove any duplicated catches at the same detector
##Pull out unique records, based on certain columns while retaining all columns
CapData<-CapData[!duplicated(CapData[,c("#Session","ID","Occassion","Detector" )]),]

##Order by bear and year
##Re-Order dataframe
CapData<-CapData[with(CapData, order(`#Session`, ID, Occassion)), ]


##Export

write.table(CapData,file=here::here("CleanData", "CapHists_SECR", "CapData", paste(df$Project[1], ".txt",sep="")), sep=",", row.names=F, quote=FALSE, col.names=TRUE)


########################################
########################################
## PREP trap data
########################################
########################################


###### Prepare Detector Layout  ##########

##Retain the columns that will be of use for Detector Layout
dfdetectcull<-df[c("SITE_NAME","Easting.DC","Northing.DC","Trap.Type","Session", "Project")]

##Add / infront of Trap.Type
dfdetectcull$Trap.Type<-paste("/",dfdetectcull$Trap.Type,sep=" ")

##Add Year_Session
dfdetectcull$Year_Session <- paste(dfdetectcull$Project,dfdetectcull$Session,sep="_")


##Prepare Usage Matrix which shows when each site was run
Use <- reshape2::dcast(dfdetectcull,SITE_NAME + Trap.Type ~ Year_Session, value.var="Project", 
                       fun.aggregate=length)



##Replace anything>1 with a 1
Use[,3:(ncol(Use))][Use[,3:(ncol(Use))]>1] <- 1

##combine usage into a single column
Usage<-unite(Use,"ch", c(3:(ncol(Use))), sep="")

## Adf Usage to data
dfdetectcull<-merge(dfdetectcull%>%dplyr::select(-Trap.Type),Usage,by="SITE_NAME",all=TRUE)

##Rename Columns
DetectLayout<-plyr::rename(dfdetectcull, c("SITE_NAME"="#Detector", 
                                           "Easting.DC"="X",
                                           "Northing.DC"="Y",
                                           "ch"="Usage",
                                           "Trap.Type" = "Trap_Type"))


##remove any duplicates
DetectLayout<-DetectLayout[!duplicated(DetectLayout[,c("#Detector")]),]


##Retain the columns that will be of use for Detector Layout
DetectLayout<-DetectLayout[c("#Detector","X","Y","Usage", "Trap_Type")]

###throw error if
if (sum(!CapData$Detector %in% DetectLayout$`#Detector`)>0){ 
  stop("failed to match some capture locations to detector sites")
}

##round coordinates to nearest meter
DetectLayout$X<- round(DetectLayout$X,0)
DetectLayout$Y <- round(DetectLayout$Y,0)
##Export 
write.table(DetectLayout,file=here::here("CleanData", "CapHists_SECR", "TrapData", paste(df$Project[1], ".txt",sep="")), sep=",", row.names=F, quote=FALSE, col.names=TRUE)


print(df$Project[1])
}








