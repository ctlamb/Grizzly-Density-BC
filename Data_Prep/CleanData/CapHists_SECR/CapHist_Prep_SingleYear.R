
#################################
## MASSAGE DATA INTO SECR FORMAT
#################################


##Load packages
library(plyr)
require(sp)
require(rgdal)
require(maps)
require(raster)
require(rgeos)
require(ggmap)
library(reshape2)
library(data.table)
library(tidyverse)

##################################################################
####LOAD CAPTURE DATA
##################################################################

temp = list.files(path=c("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/CleanData/ProjectFiles/Single_Year",
                         "/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/CleanData/ProjectFiles/Multi_Year"),pattern = "*.csv", full.names=TRUE)


all <- ldply(temp, read_csv)



########################################
## SET UP LOOP
########################################
Proj <- unique(all$Project)
for(i in 1:length(Proj)){
  
  ##pull out a single inventory
  df<-all%>%filter(Project%in%Proj[i])
  

  ##Get into BC Albers instead of UTM11
  ####Make Data Spatial
  
  ##Loop to get each zone into correct UTM's and then project to lat long
  Zones<-unique(df$Zone)
  for(j in 1:length(Zones)){
    a<-subset(df,Zone==Zones[j])
    coordinates(a) <- c("Easting.DC","Northing.DC")
    proj4string(a)<-paste("+proj=utm +zone=",Zones[j]," +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",sep="")
    # spTransform function provide transformation between projections.
    a<- spTransform(a, CRS("+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs"))%>%
      as.data.frame()%>%
      dplyr::select("Project", "SITE_NAME", "Trap.Type", "Session", "Easting.DC", "Northing.DC", "Individual", "Sex" )
    df <- a
  }
  

  
  ##standardize U and UC sex
  df$Sex <- ifelse(df$Sex%in%c("UC", "U"),NA,as.character(df$Sex ))
  
  ##Convert Individual, Sex, SITE_NAME, and Trap.Type to factors
  df$Individual <-as.character(df$Individual)
  df$Sex <-as.factor(df$Sex)
  df$Trap.Type <-as.factor(df$Trap.Type)
  
  ##Convert SITE_NAME to universal name (X_Y) then to Factor
  df$SITE_NAME <- paste(abs(round(df$Easting.DC,0)),round(df$Northing.DC,0), df$Trap.Type, sep="_")
  df$SITE_NAME <-as.factor(df$SITE_NAME)
  
  ##rename sites, as SECR doesn't like current names
  #df <- transform(df, SITE_NAME=as.factor(match(SITE_NAME, unique(SITE_NAME))))

  ##remove whitepace from project names
  df$Project <- str_replace_all(df$Project, " ", "_")


  ##address naming issue, some bears called same name between projects, but not actually same bear
    df$Individual <- ifelse(!is.na(df$Individual) & df$Individual!="",  paste(gsub('[0-9]+', '',sapply(strsplit(df$Project, "_"), function(x){
      paste(substring(x, 1, 1), collapse = "")})),  df$Individual, sep="_"),
      as.character(df$Individual))
    


########################################
## Now put into SECR format
########################################


  
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

write.table(CapData,file=paste("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/CleanData/CapHists_SECR/CapData/", df$Project[1], ".txt",sep=""), sep=",", row.names=F, quote=FALSE, col.names=TRUE)


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
write.table(DetectLayout,file=paste("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/CleanData/CapHists_SECR/TrapData/", df$Project[1], ".txt",sep=""), sep=",", row.names=F, quote=FALSE, col.names=TRUE)


print(df$Project[1])
}






