## Preparing a selection of Alterra soil profiles (www.bodemdata.nl) for SoilGrids
## by Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(sp)
library(plotKML)

SITE <- read.csv("Alterra_Soil_xy.csv")
str(SITE)
## 643 profiles
SITE$SOURCEDB = "Alterra-BODEMDATA"
SITE$SOURCEID <- paste("Alterra", SITE$Profile_id, sep="_")
SITE.s <- SITE[!duplicated(SITE$SOURCEID),]
coordinates(SITE.s) <- ~ X + Y
proj4string(SITE.s) <- "+init=epsg:28992"
plotKML(SITE.s["Land.cover"])
SITE.s <- as.data.frame(spTransform(SITE.s, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")))
SITE.s$TIMESTRR <- as.Date(paste(SITE.s$Year), format="%Y")
SITE.s <- rename(SITE.s, c("X"="LONWGS84", "Y"="LATWGS84"))
View(SITE.s)

horizons <- read.table("Alterra_Soil_data.csv", header =TRUE, na.strings = c("-99","NA"), sep=",")
str(horizons)
## Includes also the classification:
s <- summary(horizons$Soil.Classification.code, maxsum = 140)
soiltype <- data.frame(soil_name=attr(s, "names"), count=s)
write.csv(soiltype, "Soil.Classification.code_count.csv")
legNL <- read.csv("cleanup_Alterra.csv")
horizons$TAXNWRB <- join(horizons, legNL, by="Soil.Classification.code", type="left")$TAXNWRB
summary(horizons$TAXNWRB)

horizons$ORCDRC <- horizons$Organic.Matter....*10 /1.724     ## OC in permilles
#horizons$ORCDRC <- ifelse(is.na(horizons$Carbon....), horizons$ORCDRC, horizons$Carbon....*10)
summary(horizons$ORCDRC)
hist(log1p(horizons$ORCDRC))
horizons$SOURCEID <- paste("Alterra", horizons$Profile_id, sep="_")
horizons$SAMPLEID <- make.unique(paste(horizons$Profile_id, horizons$Layer, sep="_"))
horizons <- rename(horizons, c("pH.KCl"="PHIKCL", "Sand...50mu....."="SNDPPT", "Clay...2.mu....."="CLYPPT", "Silt..2.50.mu....."="SLTPPT", "CEC..mmol.kg."="CECSUM", "Bulk.density..g.per.cm3."="BLD", "Start..cm."="UHDICM", "End..cm."="LHDICM"))
summary(horizons$UHDICM)
summary(horizons$LHDICM)
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
sumTex <- rowSums(horizons[,c("SLTPPT","CLYPPT","SNDPPT")])
hist(sumTex)
## some horizons do not sum up to 100% but these do not have to be filtered out
## Convert to cmol / kg
horizons$CECSUM <- horizons$CECSUM/10

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.Alterra <- join(SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], horizons[,c("SOURCEID","TAXNWRB")], type="left", match="first")
TAXNWRB.Alterra <- TAXNWRB.Alterra[!is.na(TAXNWRB.Alterra$TAXNWRB)&!is.na(TAXNWRB.Alterra$LONWGS84)&nchar(paste(TAXNWRB.Alterra$TAXNWRB))>0,]
str(TAXNWRB.Alterra)
### 344 profiles
coordinates(TAXNWRB.Alterra) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.Alterra) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.Alterra["TAXNWRB"])
save(TAXNWRB.Alterra, file="TAXNWRB.Alterra.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.Alterra <- join(horizons[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIKCL","ORCDRC","CECSUM","BLD")], SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
SPROPS.Alterra <- SPROPS.Alterra[!is.na(SPROPS.Alterra$LONWGS84) & !is.na(SPROPS.Alterra$LATWGS84) & !is.na(SPROPS.Alterra$DEPTH),]
str(SPROPS.Alterra)
## 2617
save(SPROPS.Alterra, file="SPROPS.Alterra.rda")
plot(SPROPS.Alterra$LONWGS84, SPROPS.Alterra$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

horizons.s <- horizons[!is.na(horizons$Horizon.code)&nchar(paste(horizons$Horizon.code))>0,]
sel.r <- grep(pattern="^R", horizons.s$DESIG, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", horizons.s$DESIG, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- grep(pattern="CR", horizons.s$DESIG, ignore.case=FALSE, fixed=FALSE)
horizons.s$BDRICM <- NA
horizons.s$BDRICM[sel.r] <- horizons.s$UHDICM[sel.r]
horizons.s$BDRICM[sel.r2] <- horizons.s$LHDICM[sel.r2]
horizons.s$BDRICM[sel.r3] <- horizons.s$LHDICM[sel.r3]
bdr.d <- aggregate(horizons.s$BDRICM, list(horizons.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.Alterra <- join(SITE.s[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84")], bdr.d, type="left")
BDR.Alterra$BDRICM <- ifelse(is.infinite(BDR.Alterra$BDRICM), 250, BDR.Alterra$BDRICM)
BDR.Alterra <- BDR.Alterra[!is.na(BDR.Alterra$BDRICM),]
str(BDR.Alterra)
## 643 points --> no bedrock observed!
save(BDR.Alterra, file="BDR.Alterra.rda") 
