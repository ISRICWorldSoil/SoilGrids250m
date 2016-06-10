## Reading and writing of eSOTER DB (2213 profiles with no lab data);
## eSOTER project [http://www.esoter.net/?q=content/wp2-methodology-integrate-soil-data-legacy-and-remote-sensing-sources]
## Tom.Hengl@isric.org

library(RODBC)
library(aqp)
library(plyr)
library(sp)
library(plotKML)

# ------------------------------------------------------------
# Fetch tables
# ------------------------------------------------------------

# import to R:
cGSPD <- odbcConnect(dsn="eSOTER")
odbcGetInfo(cGSPD)
sqlTables(cGSPD)$TABLE_NAME
# get horizon table:
#SITE <- sqlFetch(cGSPD, "Profile", stringsAsFactors=FALSE, as.is=TRUE)
str(SITE)
#HORIZON <- sqlFetch(cGSPD, "RepresentativeHorizonValues", stringsAsFactors=FALSE, as.is=TRUE)  
HORIZON <- read.csv("eSOTER_RepresentativeHorizonValues.csv") 
str(HORIZON)

# ------------------------------------------------------------
# Re-format columns
# ------------------------------------------------------------

SITE <- subset(SITE, !is.na(SITE$LATI)|!is.na(SITE$LNGI))  # no point in using profiles that have no geo-reference!
## rename columns:
SITE.s <- rename(SITE, c("PRID"="SOURCEID", "LATI"="LATWGS84", "LNGI"="LONWGS84", "SAYR"="TIMESTRR", "WRBC"="TAXNWRB", "STAX"="TAXNUSDA"))
SITE.s$SOURCEDB = "eSOTER"
SITE.s$TIMESTRR <- as.Date(paste(SITE.s$TIMESTRR), format="%Y")
sel.r <- grep(pattern="leptic", SITE.s$TAXNWRB, ignore.case=TRUE, fixed=FALSE)
SITE.s$TAXGWRB <- sapply(SITE.s$TAXNWRB, function(x){strsplit(x, " ")[[1]][length(strsplit(x, " ")[[1]])]})
View(SITE.s)

HOR.s <- rename(HORIZON, c("PRID"="SOURCEID", "HBTP"="UHDICM", "HBDE"="LHDICM", "PHAQ"="PHIHOX", "SDTO"="SNDPPT", "STPC"="SLTPPT", "CLPC"="CLYPPT", "SDVC"="CRFVOL", "BULK"="BLD", "ORGC"="ORCDRC", "CECS"="CEC"))
HOR.s$ORCDRC <- HOR.s$ORCDRC*10
HOR.s$SAMPLEID <- make.unique(paste(HOR.s$SOURCEID, HOR.s$HONU, sep="_"))
HOR.s$DEPTH <- HOR.s$UHDICM + (HOR.s$LHDICM - HOR.s$UHDICM)/2
summary(HOR.s$BLD)
HOR.s$BLD <- HOR.s$BLD*1000
HOR.s$BLD[HOR.s$BLD>2600] <- NA
summary(HOR.s$ORCDRC)

## Depth to bedrock i.e. 'R' horizon (there is a column 'RockDpth' in the Profiles table but it is empty):
sel.r1 <- grep(pattern="^R", HOR.s$HODE, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", HOR.s$HODE, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", HOR.s$HODE, ignore.case=FALSE, fixed=FALSE)
HOR.s$BDRICM <- NA
## Also upper depths for R horizon entered!
HOR.s$BDRICM[sel.r1] <- HOR.s$UHDICM[sel.r1]
HOR.s$BDRICM[sel.r2] <- HOR.s$UHDICM[sel.r2]
HOR.s$BDRICM[sel.r3] <- HOR.s$UHDICM[sel.r3]
bdr.d <- aggregate(HOR.s$BDRICM, list(HOR.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.eSOTER <- join(bdr.d, SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
BDR.eSOTER$BDRICM <- ifelse(is.infinite(BDR.eSOTER$BDRICM), 250, BDR.eSOTER$BDRICM)
summary(BDR.eSOTER$BDRICM<250)
## 40 points
str(BDR.eSOTER)
save(BDR.eSOTER, file="BDR.eSOTER.rda")


# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.eSOTER <- SITE.s[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.eSOTER <- TAXNWRB.eSOTER[!is.na(TAXNWRB.eSOTER$TAXNWRB)&nchar(paste(TAXNWRB.eSOTER$TAXNWRB))>0,]
## 2213 profiles
coordinates(TAXNWRB.eSOTER) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.eSOTER) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.eSOTER["TAXNWRB"])
save(TAXNWRB.eSOTER, file="TAXNWRB.eSOTER.rda")

TAXOUSDA.eSOTER <- profs.f@site[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNUSDA")]
TAXOUSDA.eSOTER <- TAXOUSDA.eSOTER[!is.na(TAXOUSDA.eSOTER$TAXNUSDA)&nchar(paste(TAXOUSDA.eSOTER$TAXNUSDA))>0,]
## 145 profiles
coordinates(TAXOUSDA.eSOTER) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.eSOTER) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.eSOTER["TAXNUSDA"])
save(TAXOUSDA.eSOTER, file="TAXOUSDA.eSOTER.rda")

## export taxonomy data for correlation:
eSOTER_tax <- SITE.s[,c("SOURCEID", "LATWGS84", "LONWGS84", "TAXNWRB", "TAXNUSDA", "TAXGWRB", "LITH")]
summary(as.factor(eSOTER_tax$TAXNUSDA))
save(eSOTER_tax, file="eSOTER_tax.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.eSOTER <- join(HOR.s[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHOX","ORCDRC")], SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
SPROPS.eSOTER <- SPROPS.eSOTER[!is.na(SPROPS.eSOTER$LONWGS84) & !is.na(SPROPS.eSOTER$LATWGS84),]
View(SPROPS.eSOTER)
## 5899
save(SPROPS.eSOTER, file="SPROPS.eSOTER.rda")