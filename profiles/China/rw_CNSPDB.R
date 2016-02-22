## Reading and writing of China SOTER DB (1430 profiles) and National soil profile DB
## SOTER project for China [http://www.isric.org/data/soil-and-terrain-database-china]
## CSPBNU = "A China data set of soil properties for land surface modeling" doi:10.1002/jame.20026, 2013
## Profile data prepared by Wei Shangguan (https://scholar.google.com/citations?user=sWZZ984AAAAJ&hl=en)
# Tom.Hengl@isric.org

#library(RODBC)
library(aqp)
library(plyr)
library(sp)
library(plotKML)

# ------------------------------------------------------------
# Fetch tables
# ------------------------------------------------------------

## import to R:
#cGSPD <- odbcConnect(dsn="CN")
#odbcGetInfo(cGSPD)
#sqlTables(cGSPD)$TABLE_NAME
## get horizon table:
#SITE <- sqlFetch(cGSPD, "Profile", stringsAsFactors=FALSE, as.is=TRUE)
SITE <- read.csv("CHINA_SOTERv1_Profile.csv")
str(SITE)
#HORIZON <- sqlFetch(cGSPD, "RepresentativeHorizonValues", stringsAsFactors=FALSE, as.is=TRUE)  
HORIZON <- read.csv("CHINA_SOTERv1_Horizon.csv")
horizons <- read.csv("china_wosis_horizon.csv")
str(horizons)
horizons$SAMPLEID <- make.unique(paste("CN", horizons$SOURCEID, horizons$shn, sep="_"))
horizons$SOURCEID <- paste("CN", horizons$SOURCEID, sep="_")
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
horizons$BLD <- horizons$BLD*1000
## Total carbon?
horizons$ORCDRC <- horizons$ORCDRC*10

sites <- read.csv("china_wosis_sites.csv")
#sites$TIMESTRR <- as.Date(paste(sites$TIMESTR), format="%m-%d-%Y")
sites$SOURCEID <- paste("CN", sites$SOURCEID, sep="_")

# ------------------------------------------------------------
# Re-format columns
# ------------------------------------------------------------

SITE.s <- subset(SITE, !is.na(SITE$LATI)|!is.na(SITE$LNGI))  # no point in using profiles that have no geo-reference!
## rename columns:
SITE.s <- rename(SITE.s, c("PRID"="SOURCEID", "LATI"="LATWGS84", "LNGI"="LONWGS84", "CLAF"="FAO_90"))
SITE.s$SOURCEDB = "CN-SOTER"
legFAO_90 <- read.csv("cleanup_SU_SYM90.csv")
SITE.s$TAXNWRB <- join(SITE.s, legFAO_90, type="left")$Name
View(SITE.s)

HOR.s <- rename(HORIZON, c("PRID"="SOURCEID", "HBDE"="LHDICM", "PHAQ"="PHIHO5", "SDTO"="SNDPPT", "STPC"="SLTPPT", "CLPC"="CLYPPT", "SDVC"="CRFVOL", "BULK"="BLD", "TOTC"="ORCDRC", "CECS"="CEC"))
HOR.s$SOURCEID <- paste("CN", HOR.s$SOURCEID, sep="_")
## upper depth missing!
HOR.s$UHDICM <- NA
HOR.s$UHDICM[2:nrow(HOR.s)] <- HOR.s$LHDICM[1:(nrow(HOR.s)-1)]
HOR.s$UHDICM <- ifelse(HOR.s$UHDICM > HOR.s$LHDICM, 0, HOR.s$UHDICM)
summary(HOR.s$ORCDRC)
HOR.s$SNDPPT <- as.numeric(HOR.s$SNDPPT)
HOR.s$SLTPPT <- as.numeric(HOR.s$SLTPPT)
HOR.s$CLYPPT <- as.numeric(HOR.s$CLYPPT)
summary(HOR.s$SNDPPT)
summary(HOR.s$CRFVOL)
## Depth to bedrock:
#sel.r <- grep(pattern="R", HOR.s$HODE, ignore.case=FALSE, fixed=FALSE)
## only two horizons?

## hits on visual check:
## 18 obs with deeper 'upper' HOR.s
str(HOR.s[which(HOR.s$LHDICM-HOR.s$UHDICM<=0),])

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.CHINA_SOTERv1 <- SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.CHINA_SOTERv1 <- TAXNWRB.CHINA_SOTERv1[!is.na(TAXNWRB.CHINA_SOTERv1$TAXNWRB)&!is.na(TAXNWRB.CHINA_SOTERv1$LONWGS84)&nchar(paste(TAXNWRB.CHINA_SOTERv1$TAXNWRB))>0,]
str(TAXNWRB.CHINA_SOTERv1)
## 1429 profiles
coordinates(TAXNWRB.CHINA_SOTERv1) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.CHINA_SOTERv1) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.CHINA_SOTERv1["TAXNWRB"])
save(TAXNWRB.CHINA_SOTERv1, file="TAXNWRB.CHINA_SOTERv1.rda")

##visual check
## no detected issues

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.CSPBNU <- plyr::join(horizons[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHOX","BLD","ORCDRC","CECSUM")], sites[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
SPROPS.CSPBNU <- SPROPS.CSPBNU[!is.na(SPROPS.CSPBNU$LONWGS84) & !is.na(SPROPS.CSPBNU$LATWGS84) & !is.na(SPROPS.CSPBNU$DEPTH),]
#View(SPROPS.CSPBNU)
str(SPROPS.CSPBNU)
## 32208 horizons
summary(SPROPS.CSPBNU$BLD)
summary(SPROPS.CSPBNU$ORCDRC)
save(SPROPS.CSPBNU, file="SPROPS.CSPBNU.rda")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

horizons.s <- horizons[!is.na(horizons$H)&nchar(paste(horizons$H))>0,]
sel.r <- grep(pattern="^R", horizons.s$H, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", horizons.s$H, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", horizons.s$H, ignore.case=FALSE, fixed=FALSE)
horizons.s$BDRICM <- NA
horizons.s$BDRICM[sel.r] <- horizons.s$UHDICM[sel.r]
horizons.s$BDRICM[sel.r2] <- horizons.s$LHDICM[sel.r2]
horizons.s$BDRICM[sel.r3] <- horizons.s$LHDICM[sel.r3]
bdr.d <- aggregate(horizons.s$BDRICM, list(horizons.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.CSPBNU <- join(sites[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], bdr.d, type="left")
BDR.CSPBNU$BDRICM <- ifelse(is.infinite(BDR.CSPBNU$BDRICM), 250, BDR.CSPBNU$BDRICM)
BDR.CSPBNU <- BDR.CSPBNU[!is.na(BDR.CSPBNU$BDRICM),]
str(BDR.CSPBNU)
summary(BDR.CSPBNU$BDRICM<250)
## 6988 points
save(BDR.CSPBNU, file="BDR.CSPBNU.rda") 