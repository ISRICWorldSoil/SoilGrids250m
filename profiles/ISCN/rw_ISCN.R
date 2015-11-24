## ISCN soil data (http://iscn.fluxdata.org/Data/Pages/DatabaseReports.aspx)
## by: Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(plotKML)

ISCN <- read.csv("ISCNProfileData_LATEST.csv")
str(ISCN)
ISCN$SOURCEID <- paste(ISCN$profile_name)
ISCN$TIMESTRR <- as.Date(ISCN$observation_date, format="%d-%m-%Y")
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

SPROPS.ISCN <- join(as.data.frame(TAXOUSDA.ISCN)[,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB")], HORIZON[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CRFVOL","SNDPPT","CLYPPT","BLD","SLTPPT","PHIHOX","ORCDRC")], type="left")
SPROPS.ISCN <- SPROPS.ISCN[!is.na(SPROPS.ISCN$DEPTH),]
str(SPROPS.ISCN)
summary(SPROPS.ISCN$ORCDRC) ## mean=6.9% !!
summary(SPROPS.ISCN$BLD) ## mean=1381
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
