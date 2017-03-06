## ISCN soil data (http://iscn.fluxdata.org/Data/Pages/DatabaseReports.aspx)
## by: Tom.Hengl@isric.org

load(".RData")
library(aqp)
library(plyr)
library(stringr)
library(sp)
library(plotKML)

ISCN <- read.csv("ISCNProfileData_LATEST.csv")
str(ISCN)
ISCN$SOURCEID <- paste(ISCN$profile_name)
ISCN$TIMESTRR <- as.Date(ISCN$observation_date, format="%m/%d/%Y")
summary(ISCN$TIMESTRR)
summary(ISCN$contributorOrganization)
## soil taxonomy unstandardized:
ISCN$TAXNUSDA <- ifelse(is.na(ISCN$soil_taxon), paste(ISCN$soil_series), paste(ISCN$soil_taxon))
ISCN$SOURCEDB <- gsub(" ", "_", ISCN$contributorOrganization)

## Mask out NRCS data because the original source is more reliable
selP <- !ISCN$contributorOrganization=="NRCS Jan/2011"
TAXOUSDA.ISCN <- ISCN[selP,c("SOURCEID","lat..dec..deg.","long..dec..deg.","datum","TIMESTRR","TAXNUSDA","SOURCEDB")]
summary(TAXOUSDA.ISCN$datum)
## conversion of coordinates is required:
TAXOUSDA.ISCN$LONWGS84 <- TAXOUSDA.ISCN$long..dec..deg. 
TAXOUSDA.ISCN$LATWGS84 <- TAXOUSDA.ISCN$lat..dec..deg. 
TAXOUSDA.ISCN <- TAXOUSDA.ISCN[TAXOUSDA.ISCN$LONWGS84>-180&TAXOUSDA.ISCN$LONWGS84<180&TAXOUSDA.ISCN$LATWGS84>-90&TAXOUSDA.ISCN$LATWGS84<90,]
str(TAXOUSDA.ISCN)
## 3401 profiles

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

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

ISCN.hor <- read.csv("ISCNLayerData_LATEST.csv")
summary(ISCN.hor$oc....)
HORIZON <- ISCN.hor[!ISCN.hor$contributorOrganization=="NRCS Jan/2011",c("profile_name", "layer_name", "layer_top..cm.", "layer_bot..cm.", "hzn_desgn", "oc....", "bd_samp..g.cm.3.", "ph_h2o", "sand_tot_psa....", "silt_tot_psa....", "clay_tot_psa....", "wpg2....")]
HORIZON$SOURCEID <- paste(HORIZON$profile_name)
HORIZON$SAMPLEID <- make.unique(paste(HORIZON$profile_name, HORIZON$layer_name, sep="_"))
HORIZON$SOURCEDB = "ISCN"
HORIZON <- rename(HORIZON, c("layer_top..cm."="UHDICM", "layer_bot..cm."="LHDICM", "oc...."="ORCDRC", "bd_samp..g.cm.3."="BLD", "ph_h2o"="PHIHOX", "sand_tot_psa...."="SNDPPT", "silt_tot_psa...."="SLTPPT", "clay_tot_psa...."="CLYPPT", "wpg2...."="CRFVOL"))
HORIZON$DEPTH <- HORIZON$UHDICM + (HORIZON$LHDICM - HORIZON$UHDICM)/2
HORIZON$BLD <- HORIZON$BLD * 1000
HORIZON$BLD <- ifelse(HORIZON$BLD<100, NA, HORIZON$BLD)
summary(HORIZON$BLD)
HORIZON$ORCDRC <- HORIZON$ORCDRC * 10
summary(HORIZON$ORCDRC)
HORIZON$PHIHOX <- ifelse(HORIZON$PHIHOX>11, HORIZON$PHIHOX/10, HORIZON$PHIHOX)
HORIZON$PHIHOX <- ifelse(HORIZON$PHIHOX<2.5, NA, HORIZON$PHIHOX)
summary(HORIZON$PHIHOX)

SPROPS.ISCN <- join(as.data.frame(TAXOUSDA.ISCN)[,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB","TIMESTRR")], HORIZON[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CRFVOL","SNDPPT","CLYPPT","BLD","SLTPPT","PHIHOX","ORCDRC")], type="left")
SPROPS.ISCN <- SPROPS.ISCN[!is.na(SPROPS.ISCN$DEPTH),]
str(SPROPS.ISCN)
summary(SPROPS.ISCN$ORCDRC) ## mean=6.9% !!
summary(SPROPS.ISCN$BLD) ## mean=705
## 9450
save(SPROPS.ISCN, file="SPROPS.ISCN.rda")
plot(SPROPS.ISCN$LONWGS84, SPROPS.ISCN$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

HORIZON.s <- HORIZON[!is.na(HORIZON$hzn_desgn)&nchar(paste(HORIZON$hzn_desgn))>0,]
summary(HORIZON.s$hzn_desgn)
sel.r <- grep(pattern="^R", as.character(HORIZON.s$hzn_desgn), ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", as.character(HORIZON.s$hzn_desgn), ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", as.character(HORIZON.s$hzn_desgn), ignore.case=FALSE, fixed=FALSE)
HORIZON.s$BDRICM <- NA
HORIZON.s$BDRICM[sel.r] <- HORIZON.s$UHDICM[sel.r]
HORIZON.s$BDRICM[sel.r2] <- HORIZON.s$LHDICM[sel.r2]
HORIZON.s$BDRICM[sel.r3] <- HORIZON.s$LHDICM[sel.r3]
bdr.d <- aggregate(HORIZON.s$BDRICM, list(HORIZON.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.ISCN <- join(as.data.frame(TAXOUSDA.ISCN)[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], bdr.d, type="left")
BDR.ISCN$BDRICM <- ifelse(is.infinite(BDR.ISCN$BDRICM), 250, BDR.ISCN$BDRICM)
BDR.ISCN <- BDR.ISCN[!is.na(BDR.ISCN$BDRICM),]
str(BDR.ISCN)
summary(BDR.ISCN$BDRICM<250)
## only 7 points!
save(BDR.ISCN, file="BDR.ISCN.rda") 

# ------------------------------------------------------------
# Soil organic carbon stock (kg / m2)
# ------------------------------------------------------------

SOCS.ISCN = read.csv("ISCNData_12_Dec_2016_Profile.csv")
SOCS.ISCN$SOURCEID <- paste(SOCS.ISCN$X.profile_name.)
SOCS.ISCN$TIMESTRR <- as.Date(paste(SOCS.ISCN$X.observation_date..YYYY.MM.DD..), format="%Y-%m-%d")
summary(SOCS.ISCN$TIMESTRR)
SOCS.ISCN$SOURCEDB = as.factor(SOCS.ISCN$X.dataset_name_sub.)
summary(SOCS.ISCN$SOURCEDB)
## Mask out NRCS data because the original source is more reliable
selP2 <- !SOCS.ISCN$SOURCEDB=="NRCS Sept/2014"
SOCS.ISCN <- SOCS.ISCN[selP2,c("SOURCEID","X.lat..dec..deg..","X.long..dec..deg..","X.datum..datum..","TIMESTRR","SOURCEDB","X.soc..g.cm.2..","X.soc_method.","X.soc_depth..cm..")]
SOCS.ISCN = SOCS.ISCN[!is.na(SOCS.ISCN$X.lat..dec..deg..),]
SOCS.ISCN$YEAR = as.numeric(format(SOCS.ISCN$TIMESTRR,'%Y'))
## 19,016 points
summary(SOCS.ISCN$X.datum..datum..)
SOCS.ISCN$LONWGS84 = SOCS.ISCN$X.long..dec..deg..
SOCS.ISCN$LATWGS84 = SOCS.ISCN$X.lat..dec..deg..
for(j in c("NAD83","NAD27")){
  sel <- which(SOCS.ISCN$X.datum..datum..==j)
  xy <- SOCS.ISCN[sel,c("X.long..dec..deg..","X.lat..dec..deg..")]
  coordinates(xy) <- ~ X.long..dec..deg.. + X.lat..dec..deg..
  if(j=="NAD 83"|j=="NAD83"|j=="NAD83?"){
    proj4string(xy) <- CRS("+proj=longlat +datum=NAD83")
  }
  if(j=="NAD27"){
    proj4string(xy) <- CRS("+proj=longlat +datum=NAD27")
  }
  xy <- spTransform(xy, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
  SOCS.ISCN[sel,"LONWGS84"] = xy@coords[,1]
  SOCS.ISCN[sel,"LATWGS84"] = xy@coords[,2]
}

SOCS.ISCN$dSOCS_100cm = NA
SOCS.ISCN$dSOCS_200cm = NA
SOCS.ISCN[SOCS.ISCN$X.soc_depth..cm..==100,"dSOCS_100cm"] = SOCS.ISCN[SOCS.ISCN$X.soc_depth..cm..==100,"X.soc..g.cm.2.."]*10
SOCS.ISCN[SOCS.ISCN$X.soc_depth..cm..==200,"dSOCS_200cm"] = SOCS.ISCN[SOCS.ISCN$X.soc_depth..cm..==200,"X.soc..g.cm.2.."]*10
SOCS.ISCN = SOCS.ISCN[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","dSOCS_100cm","dSOCS_200cm","YEAR")]
summary(SOCS.ISCN$dSOCS_100cm)
summary(SOCS.ISCN$dSOCS_200cm)
coordinates(SOCS.ISCN) <- ~ LONWGS84 + LATWGS84
proj4string(SOCS.ISCN) = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
plot(SOCS.ISCN)
save(SOCS.ISCN, file="SOCS.ISCN.rda")

