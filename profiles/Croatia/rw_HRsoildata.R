# title         : rw_HRsoildata.R
# purpose       : Reading and writing of the Croatian soil profile data;
# reference     : The Croatian national soil profile db (2199 profiles);
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Aug 2012.
# inputs        : "HRSGDB_profiles.txt" all in one table
# outputs       : Horizons and sites data frames;
# remarks 1     : the table data have been published in a book -- Martinovic J, Vrankovic A (1997). Baza tala Republike Hrvatske. I-III. Ministarstvo zaštite okoliša i prostornog planiranja, Zagreb. p. 365., but the data is not publicly available;

#download.file("http://globalsoilmap.net/data/HRSGDB_profiles.zip", destfile=paste(getwd(),"HRSGDB_profiles.zip",sep="/"))
#unzip("HRSGDB_profiles.zip")

library(aqp)
library(plyr)
library(rgdal)
library(GSIF)

## Import the soil profile data HRSDGB:
HRSGDB <- read.delim("HRSGDB_profiles.txt", stringsAsFactors = FALSE)
str(HRSGDB)
HRSGDB$corrX <- ifelse(!is.na(HRSGDB$corrX), HRSGDB$corrX, HRSGDB$Cro16.30_X)
HRSGDB$corrY <- ifelse(!is.na(HRSGDB$corrY), HRSGDB$corrY, HRSGDB$Cro16.30_Y)
HRSGDB <- HRSGDB[!is.na(HRSGDB$corrX)&!is.na(HRSGDB$corrY),]
coordinates(HRSGDB) <- ~ corrX + corrY
HRSGDB$Cro16.30_X <- NULL
HRSGDB$Cro16.30_Y <- NULL
proj4string(HRSGDB) <- CRS("+proj=tmerc +lat_0=0 +lon_0=16.5 +k=0.9999 +x_0=2500000 +y_0=0 +ellps=bessel +towgs84=550.499,164.116,475.142,5.80967,2.07902,-11.62386,0.99999445824 +units=m")
HRSGDB.ll <- spTransform(HRSGDB, CRS("+proj=longlat +datum=WGS84"))
site <- data.frame(HRSGDB.ll)

# rename columns:
names(site)[which(names(site)=="corrX")] <- "LONWGS84" # longitude
names(site)[which(names(site)=="corrY")] <- "LATWGS84" # latitude
site$SOURCEID <- as.factor(paste("HR", site$PROF_ID, sep="_"))
site$TIMESTRR <- as.Date(as.character(site$UZORAK), format="%Y")
site$SOURCEDB <- "HRSPDB"
FAOleg <- read.csv("../Full_data_FAO74_US_CPSS.csv")
FAOleg <- FAOleg[!is.na(FAOleg$SoilUnitD90)&nchar(paste(FAOleg$SoilUnitD90))>0,]
FAOleg$FAO_NAME <- FAOleg$FAO74
site$TAXNWRB <- join(site, FAOleg, type="left", match = "first")$SoilUnitD90

## sort horizons:
horizon <- getHorizons(site, idcol="SOURCEID", sel=c("OZN", "GOR", "DON", "MKP", "PH1", "PH2", "MSP", "MP", "MG", "HUM"), pattern=paste0("HOR", 1:9,"_"), reverse=TRUE)
# TH: Silt is in 0.02-0.002 mm and not in the 0.05-0.002 system!!
horizon$UHDICM <- as.numeric(horizon$GOR)
horizon$LHDICM <- as.numeric(horizon$DON)
horizon$PHIHOX <- as.numeric(horizon$PH1) 
horizon$PHIKCL <- as.numeric(horizon$PH2)
horizon$SNDPPT <- as.numeric(horizon$MSP)*0.8 + as.numeric(horizon$MKP)
horizon$SLTPPT <- as.numeric(horizon$MP) + as.numeric(horizon$MSP)*0.2
horizon$CLYPPT <- as.numeric(horizon$MG)
horizon$ORCDRC <- as.numeric(horizon$HUM)*10/1.724
horizon$CRFVOL <- join(horizon, site[,c("SOURCEID","STIJENA")], type="left")$STIJENA
horizon$LHDICM <- ifelse(is.na(horizon$LHDICM), horizon$UHDICM+50, horizon$LHDICM)
horizon$DEPTH <- horizon$UHDICM + (horizon$LHDICM - horizon$UHDICM)/2
horizon <- horizon[!is.na(horizon$DEPTH),]
## reduce coarse fragments based for the top horizon:
horizon$CRFVOL[horizon$UHDICM<30] <- horizon$CRFVOL[horizon$UHDICM<30]*.3
horizon$SAMPLEID <- make.unique(paste(horizon$SOURCEID, horizon$OZN, sep="_"))

# filter strange values:
summary(horizon$PHIHOX)
horizon$PHIHOX <- ifelse(horizon$PHIHOX < 2 | horizon$PHIHOX > 10, NA, horizon$PHIHOX)
horizon$PHIKCL <- ifelse(horizon$PHIKCL < 2 | horizon$PHIKCL > 10, NA, horizon$PHIKCL)
horizon$SNDPPT[horizon$SAMPLEID=="HR_id 805_Amo"] <- horizon$SNDPPT[horizon$SAMPLEID=="HR_id 805_Amo"]/10
summary(horizon$ORCDRC)
## Correct texture fractions:
sumTex <- rowSums(horizon[,c("SLTPPT","CLYPPT","SNDPPT")])
horizon$SNDPPT <- horizon$SNDPPT / ((sumTex - horizon$CLYPPT) /(100 - horizon$CLYPPT))
horizon$SLTPPT <- horizon$SLTPPT / ((sumTex - horizon$CLYPPT) /(100 - horizon$CLYPPT))

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.HRSPDB <- site[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.HRSPDB <- TAXNWRB.HRSPDB[!is.na(TAXNWRB.HRSPDB$TAXNWRB)&!is.na(TAXNWRB.HRSPDB$LONWGS84)&nchar(paste(TAXNWRB.HRSPDB$TAXNWRB))>0,]
str(TAXNWRB.HRSPDB)
## 1866 profiles
coordinates(TAXNWRB.HRSPDB) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.HRSPDB) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.HRSPDB["TAXNWRB"])
save(TAXNWRB.HRSPDB, file="TAXNWRB.HRSPDB.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.HRSPDB <- join(horizon[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","CRFVOL","PHIHOX","ORCDRC")], site[,c("SOURCEID","LONWGS84","LATWGS84")], type="left")
SPROPS.HRSPDB <- SPROPS.HRSPDB[!is.na(SPROPS.HRSPDB$LONWGS84) & !is.na(SPROPS.HRSPDB$LATWGS84),]
View(SPROPS.HRSPDB)
## 5899
save(SPROPS.HRSPDB, file="SPROPS.HRSPDB.rda")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

levels(site$TAXNWRB)
horizons.s <- horizon[!is.na(horizon$OZN)&nchar(paste(horizon$OZN))>0,]
sel.r <- grep(pattern="^R", horizons.s$OZN, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", horizons.s$OZN, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", horizons.s$OZN, ignore.case=FALSE, fixed=FALSE)
sel.r4 <- unique(c(which(horizons.s$SOURCEID %in% site$SOURCEID[grep(pattern="LEPT", ignore.case=TRUE, paste(site$TAXNWRB))]), which(horizons.s$SOURCEID %in% site$SOURCEID[grep(pattern="LITH", ignore.case=TRUE, paste(site$TAXNWRB))])))
horizons.s$BDRICM <- NA
horizons.s$BDRICM[sel.r] <- horizons.s$UHDICM[sel.r]
horizons.s$BDRICM[sel.r2] <- horizons.s$LHDICM[sel.r2]
horizons.s$BDRICM[sel.r3] <- horizons.s$LHDICM[sel.r3]
horizons.s$BDRICM[sel.r4] <- horizons.s$LHDICM[sel.r4]
bdr.d <- aggregate(horizons.s$BDRICM, list(horizons.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.HRSPDB <- join(site[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], bdr.d, type="left")
BDR.HRSPDB$BDRICM <- ifelse(is.infinite(BDR.HRSPDB$BDRICM), 250, BDR.HRSPDB$BDRICM)
BDR.HRSPDB <- BDR.HRSPDB[!is.na(BDR.HRSPDB$BDRICM),]
str(BDR.HRSPDB)
summary(BDR.HRSPDB$BDRICM<250)
## 18 points
save(BDR.HRSPDB, file="BDR.HRSPDB.rda")