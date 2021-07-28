library(plyr)
library(EasyMARK)
##Load Spatial packages
require(sp)
require(rgdal)
require(maps)
require(raster)
require(rgeos)
require(maptools)
require(RgoogleMaps)
require(ggmap)
##Here is a list of general formatting requirements:
##   2)All lines in the encounter histories file must be the same length
##   3)Comments can be added to the file using this format-    /*comment */
##   5)Any covariates added to the input file should be scaled between 0-1


##Load Data

setwd("/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/RawData/From_BMcLellan")
Raw <- read.csv("2016_09_19_CoastData.csv", header=TRUE, sep=",")



###FIX ERRORS
Raw[Raw$NEW.SITE.ID=="HT Cadwalider 2010" & Raw$YR==2010 & Raw$Session==3, "Session"] <- 2 ##only one session 3, change to session 2
Raw[Raw$NEW.SITE.ID=="HT Cadwalider 2010" & Raw$YR==2010 , "Session"]

Raw[Raw$NEW.SITE.ID=="RT L McGill Tr Big Trees" & Raw$YR==2010 & Raw$Session==3, "Session"] <- 2 ##only one session 3, change to session 2
Raw[Raw$NEW.SITE.ID=="RT L McGill Tr Big Trees" & Raw$YR==2010 , "Session"]


Raw[Raw$NEW.SITE.ID=="HT Anderson" & Raw$YR==2013 & Raw$Session==5, "Session"] <- 4 ##only one session 5, change to session 4
Raw[Raw$NEW.SITE.ID=="HT Anderson" & Raw$YR==2013 , "Session"]


Raw[Raw$NEW.SITE.ID=="RT Texas 2014" & Raw$YR==2015 & Raw$Session==6, "Session"] <- 5 ##only one session 6, change to session 5
Raw[Raw$NEW.SITE.ID=="RT Texas 2014" & Raw$YR==2015 , "Session"]


Raw[Raw$NEW.SITE.ID %in% c("RT Duf L", "RT L Phelix")  & Raw$YR==2011 & Raw$Session==1, "Session"] <- 2 ##only two session 1, change to session 2
Raw[Raw$NEW.SITE.ID %in% c("RT Duf L", "RT L Phelix")  & Raw$YR==2011, "Session"]


Raw[Raw$Site.Name %in% c("RT Cad Rd 18K", "RT Anderson1","RT Seaton Rd")  & Raw$YR==2013 & Raw$Session==5, "Session"] <- 4 ##only three session 5, change to session 4
Raw[Raw$Site.Name == "RT Anderson"  & Raw$YR==2013, "Session"]


Raw[Raw$HTorRT=="HT" & Raw$YR==2014 & Raw$Session==5, "Session"] <- 4  ##many sites run, no bears caught, join w/ session 4
Raw[Raw$HTorRT=="HT" & Raw$YR==2014, "Session"]



###change names to match w/ michelle's individual covariates
Raw$NEWIDS <-as.character(Raw$NEWIDS)
Raw[Raw$NEWIDS=="Rock " , "NEWIDS"] <-"Rock"
Raw[Raw$NEWIDS=="Carl " , "NEWIDS"] <-"Carl"
Raw[Raw$NEWIDS=="F132_M221_2014_Boy" , "NEWIDS"] <-"F132/M221_2014_Boy"
Raw[Raw$NEWIDS=="Flats  " , "NEWIDS"] <-"Flats"
Raw[Raw$NEWIDS=="Georgia_Jim_2015_girl" , "NEWIDS"] <-"Georgia/Jims_2015daugher"
Raw[Raw$NEWIDS=="Haily_ Hurly" , "NEWIDS"] <-"Haily of the Hurly"
Raw[Raw$NEWIDS=="Hazel_Hawthorne" , "NEWIDS"] <-"Hazel of Hawthorne"
Raw[Raw$NEWIDS=="Loops" , "NEWIDS"] <-"Nancy_Huck_girl_Loops"
Raw[Raw$NEWIDS=="Molly_2014_girl" , "NEWIDS"] <-"Molly_2014_cub"


##change NONE in NEWIDS to NA

Raw[Raw$NEWIDS=="NONE","NEWIDS"] <- NA



###Remove capture records for now, retain RT's and HT's
summary(Raw$HTorRT)

df<-subset(Raw,HTorRT %in% c("HT","RT"))





###Remove years of data that Bruce says are just Rochertta's "Monitoring"
##Good years are:
# So, SN, 2005, 2010, 2011, 2012, 2013, 2014, 2015
# 
# and SC 2006, 2010, 2011, 2012, 2013, 2014, 2015

df<-rbind(subset(df,GBPU=="SN"&YR %in% c("2010", "2011", "2012", "2013", "2014", "2015" )),
          subset(df,GBPU=="SC"&YR %in% c("2010", "2011", "2012", "2013", "2014", "2015" ))  
          )


##Remove any missing data
df<-subset(df,!is.na(YR))

##remove session=0 as these are samples that could be from last year as 
##sites weren't burned yet.   
df<-subset(df,Session!=0)



##do some column renaming
FINAL<-
  df%>%
  mutate(Project=paste("Stein_Nahat", YR, sep="_"),
         Zone=10)%>%
  rename(Year=YR,
         SITE_NAME=Site.Name,
         Trap.Type=HTorRT,
         Easting.DC=Easting,
         Northing.DC=Northing,
         Individual=ID
         )


##Reorder
FINAL <- FINAL[, c("Project","Year", "SITE_NAME", "Trap.Type",  "Session",     "Easting.DC",  "Northing.DC", "Zone", "Individual"  , "Sex"   )]
FINAL<-FINAL[with(FINAL, order(Year, SITE_NAME, Session)), ]


##remove multiple records from same check
FINAL <- FINAL%>%
  distinct()



##EXPORT
write.csv(FINAL, "/Users/clayton.lamb/Google Drive/Documents/University/U_A/Analyses/BC_Wide_PhD/Prov_Grizz_density_oSCR/Data_Prep/CleanData/ProjectFiles/Multi_Year/Stein_Nahat_2010_2015.csv", row.names = FALSE, na="")



