## ISCN soil data (http://iscn.fluxdata.org/Data/Pages/DatabaseReports.aspx)
## by: Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(plotKML)

ISCN <- read.csv("ISCNProfileData_LATEST.csv")
str(ISCN)
ISCN$SOURCEID <- make.unique(paste(ISCN$profile_name))
ISCN$TIMESTRR <- as.Date(ISCN$observation_date, format="%d-%m-%Y")
summary(ISCN$contributorOrganization)
## soil taxonomy unstandardized:
ISCN$TAXNUSDA <- ifelse(is.na(ISCN$soil_taxon), paste(ISCN$soil_series), paste(ISCN$soil_taxon))
ISCN$SOURCEDB <- gsub(" ", "_", ISCN$contributorOrganization)

TAXOUSDA.ISCN <- ISCN[!ISCN$contributorOrganization=="NRCS Jan/2011",c("SOURCEID","lat..dec..deg.","long..dec..deg.","datum","TIMESTRR","TAXNUSDA","SOURCEDB")]
summary(TAXOUSDA.ISCN$datum)
## conversion of coordinates is required:
TAXOUSDA.ISCN$LONWGS84 <- TAXOUSDA.ISCN$long..dec..deg. 
TAXOUSDA.ISCN$LATWGS84 <- TAXOUSDA.ISCN$lat..dec..deg. 
TAXOUSDA.ISCN <- TAXOUSDA.ISCN[TAXOUSDA.ISCN$LONWGS84>-180&TAXOUSDA.ISCN$LONWGS84<180&TAXOUSDA.ISCN$LATWGS84>-90&TAXOUSDA.ISCN$LATWGS84<90,]
str(TAXOUSDA.ISCN)

for(j in c("NAD 83","NAD83","NAD83?","NAD27")){
  sel <- which(TAXOUSDA.ISCN$datum==j)
  xy <- TAXOUSDA.ISCN[sel,c("lat..dec..deg.","long..dec..deg.")]
  coordinates(xy) <- ~ long..dec..deg. + lat..dec..deg.
  if(j=="NAD 83"|j=="NAD83"|j=="NAD83?"){
    proj4string(xy) <- CRS("+proj=longlat +datum=NAD83")
  }
  if(j=="NAD27"){
    proj4string(xy) <- CRS("+proj=longlat +datum=NAD27")
  }
  xy <- spTransform(xy, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  TAXOUSDA.ISCN[sel,"LONWGS84"] = xy@coords[,1]
  TAXOUSDA.ISCN[sel,"LATWGS84"] = xy@coords[,2]
}

TAXOUSDA.ISCN <- TAXOUSDA.ISCN[!is.na(TAXOUSDA.ISCN$TAXNUSDA)&nchar(paste(TAXOUSDA.ISCN$TAXNUSDA))>0,c("SOURCEID","SOURCEDB","TIMESTRR","TAXNUSDA","LONWGS84","LATWGS84")]
#View(TAXOUSDA.ISCN)
## 1549 profiles with coordinates and taxonomy
## convert to SPDF:
coordinates(TAXOUSDA.ISCN) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.ISCN) <- "+proj=longlat +datum=WGS84"
#plotKML(TAXOUSDA.ISCN["TAXNUSDA"])
save(TAXOUSDA.ISCN, file="TAXOUSDA.ISCN.rda")


