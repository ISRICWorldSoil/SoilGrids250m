# Reading and writing of the Brazilian national soil profile data (5781 profiles) RADAMBRASIL project;
# by Tom.Hengl@isric.org

library(rgdal)
library(plotKML)
library(plyr)

## Conjunto de Datos de Perfiles de Suelos, Escala 1:250 000 Serie II (Continuo Nacional)
edaf <- readOGR("profile_brazil.shp", "profile_brazil")
str(edaf, max.level=2)
proj4string(edaf) = "+proj=longlat +datum=WGS84"
## 11,232 points!
edaf$SOURCEID <- paste(edaf$Source, edaf$OrgProfID, sep="_")
write.csv(as.data.frame(edaf), "profile_brazil.csv")
SITE <- data.frame(edaf[!duplicated(edaf$SOURCEID),c("SOURCEID","PubYear","SoilClass0")])
s <- summary(SITE$SoilClass0, maxsum=210)
soiltype <- data.frame(SoilClass0=attr(s, "names"), count=s)
write.csv(soiltype, "soiltype_count.csv")
SITE$SOURCEDB = "RadamBrasil"
legFAO_90 <- read.csv("cleanup_RadamBrasil.csv", fileEncoding="UTF-8")
SITE$TAXNWRB <- join(SITE, legFAO_90, type="left")$FAO_1990_subgroup
SITE$TAXOUSDA <- join(SITE, legFAO_90, type="left")$USDA_suborder
SITE$TIMESTRR <- as.Date(paste(SITE$PubYear), format="%Y")
SITE <- rename(SITE, c("coords.x1"="LONWGS84", "coords.x2"="LATWGS84"))
View(SITE)

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.SolosBR <- SITE[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB","TIMESTRR")]
TAXNWRB.SolosBR <- TAXNWRB.SolosBR[!is.na(TAXNWRB.SolosBR$TAXNWRB)&!is.na(TAXNWRB.SolosBR$LONWGS84)&nchar(paste(TAXNWRB.SolosBR$TAXNWRB))>0,]
str(TAXNWRB.SolosBR)
## 5,524 profiles
coordinates(TAXNWRB.SolosBR) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.SolosBR) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.SolosBR["TAXNWRB"])
save(TAXNWRB.SolosBR, file="TAXNWRB.SolosBR.rda")

TAXOUSDA.SolosBR <- SITE[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXOUSDA","TIMESTRR")]
TAXOUSDA.SolosBR <- TAXOUSDA.SolosBR[!is.na(TAXOUSDA.SolosBR$TAXOUSDA)&!is.na(TAXOUSDA.SolosBR$LONWGS84)&nchar(paste(TAXOUSDA.SolosBR$TAXOUSDA))>0,]
str(TAXOUSDA.SolosBR)
## 5,524 profiles
coordinates(TAXOUSDA.SolosBR) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.SolosBR) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.SolosBR["TAXOUSDA"])
save(TAXOUSDA.SolosBR, file="TAXOUSDA.SolosBR.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

horizons <- as.data.frame(edaf)[,c("SOURCEID","HzSimb","HzDeIn","HzDeFn","Sand","Silt","Clay","CG","pH_H2O","C","CEC_pH7","coords.x1","coords.x2")]
horizons$SAMPLEID <- make.unique(paste(horizons$SOURCEID, horizons$HzSimb, sep="_"))
horizons <- rename(horizons, c("HzDeIn"="UHDICM","HzDeFn"="LHDICM","Sand"="SNDPPT","Silt"="SLTPPT","Clay"="CLYPPT","CG"="CRFVOL","pH_H2O"="PHIHOX","C"="ORCDRC","CEC_pH7"="CECSUM","coords.x1"="LONWGS84","coords.x2"="LATWGS84"))
horizons$ORCDRC <- horizons$ORCDRC*10
summary(horizons$ORCDRC)
## filter out all zeros!
horizons[which(horizons$PHIHOX==2.5),] ## ??
horizons$PHIHOX[horizons$PHIHOX<2|horizons$PHIHOX>11] <- NA
summary(horizons$PHIHOX)
summary(horizons$SNDPPT)
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
horizons$SOURCEDB = "RadamBrasil"
nrow(horizons)
## 11232

SPROPS.SolosBR <- horizons[!is.na(horizons$DEPTH),c("SOURCEID","SAMPLEID","SOURCEDB","UHDICM","LHDICM","DEPTH","CLYPPT","CRFVOL","SNDPPT","SLTPPT","PHIHOX","ORCDRC","CECSUM","LONWGS84","LATWGS84")]
str(SPROPS.SolosBR)
## 11,232
save(SPROPS.SolosBR, file="SPROPS.SolosBR.rda")
plot(SPROPS.SolosBR$LONWGS84, SPROPS.SolosBR$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

horizons.s <- horizons[!is.na(horizons$HzSimb)&nchar(paste(horizons$HzSimb))>0,]
sel.r <- grep(pattern="^R", horizons.s$HzSimb, ignore.case=FALSE, fixed=FALSE)
sel.r2 <- grep(pattern="*/R", horizons.s$HzSimb, ignore.case=FALSE, fixed=FALSE)
sel.r3 <- c(grep(pattern="IIR", horizons.s$HzSimb, ignore.case=FALSE, fixed=FALSE), grep(pattern="CR", horizons.s$HzSimb, ignore.case=FALSE, fixed=FALSE), grep(pattern="AR", horizons.s$HzSimb, ignore.case=FALSE, fixed=FALSE))
## 'Lithic' soils:
sel.r4 <- which(horizons.s$SOURCEID %in% edaf$SOURCEID[grep(pattern="lit", ignore.case=TRUE, paste(edaf$SoilClass))])
horizons.s$BDRICM <- NA
horizons.s$BDRICM[sel.r] <- horizons.s$UHDICM[sel.r]
horizons.s$BDRICM[sel.r2] <- horizons.s$LHDICM[sel.r2]
horizons.s$BDRICM[sel.r3] <- horizons.s$LHDICM[sel.r3]
horizons.s$BDRICM[sel.r4] <- horizons.s$LHDICM[sel.r4]
bdr.d <- aggregate(horizons.s$BDRICM, list(horizons.s$SOURCEID), max, na.rm=TRUE)
names(bdr.d) <- c("SOURCEID", "BDRICM")
BDR.SolosBR <- join(SITE[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84")], bdr.d, type="left")
BDR.SolosBR$BDRICM <- ifelse(is.infinite(BDR.SolosBR$BDRICM), 250, BDR.SolosBR$BDRICM)
BDR.SolosBR <- BDR.SolosBR[!is.na(BDR.SolosBR$BDRICM),]
str(BDR.SolosBR)
summary(BDR.SolosBR$BDRICM<250)
## 287 points with R horizon or 'Litolic' classification
save(BDR.SolosBR, file="BDR.SolosBR.rda") 