# purpose       : Reading and writing of OFRA legacy data
# reference     : OFRA data available for download [http://ec2-54-93-187-255.eu-central-1.compute.amazonaws.com/]
# producer      : by Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(sp)
library(plotKML)

SITE <- read.csv("OFRA_legacy_data.csv")
str(SITE)
## 7954 layers
s <- summary(SITE$soiltype)
soiltype <- data.frame(soil_name=attr(s, "names"), count=s)
write.csv(soiltype, "soiltype_count.csv")
SITE$SOURCEID = paste(SITE$country, SITE$distric, SITE$center, SITE$cyear, SITE$longitudes, SITE$latitudes, sep="_")
SITE$SOURCEDB = "OFRA"
SITE.s <- SITE[!duplicated(SITE$SOURCEID),c("SOURCEID","SOURCEDB","longitudes","latitudes","soiltype","cyear")]
legFAO_90 <- read.csv("cleanup_OFRA.csv")
SITE.s$TAXNWRB <- join(SITE.s, legFAO_90, type="left")$Name
SITE.s$TIMESTRR <- as.Date(paste(SITE.s$cyear), format="%Y")
SITE.s <- rename(SITE.s, c("longitudes"="LONWGS84", "latitudes"="LATWGS84"))
View(SITE.s)

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.OFRA <- SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.OFRA <- TAXNWRB.OFRA[!is.na(TAXNWRB.OFRA$TAXNWRB)&!is.na(TAXNWRB.OFRA$LONWGS84)&nchar(paste(TAXNWRB.OFRA$TAXNWRB))>0,]
str(TAXNWRB.OFRA)
## 101 profiles
coordinates(TAXNWRB.OFRA) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.OFRA) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.OFRA["TAXNWRB"])
save(TAXNWRB.OFRA, file="TAXNWRB.OFRA.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

horizons <- SITE[,c("SOURCEID","som","bulkdensity","cec","soilwaterph","sandperc","siltperc","clayperc")]
horizons$SAMPLEID <- make.unique(paste(horizons$SOURCEID, "1", sep="_"))
horizons <- rename(horizons, c("sandperc"="SNDPPT","siltperc"="SLTPPT","clayperc"="CLYPPT","bulkdensity"="BLD","soilwaterph"="PHIHOX","som"="ORCDRC","cec"="CECSUM"))
horizons$ORCDRC <- horizons$ORCDRC*10/1.724
summary(horizons$ORCDRC)
horizons$BLD <- horizons$BLD * 1000
summary(horizons$BLD)
horizons$SLTPPT <- ifelse(is.na(horizons$SLTPPT), 100-(horizons$CLYPPT+horizons$SNDPPT), horizons$SLTPPT)
horizons$SNDPPT <- ifelse(is.na(horizons$SNDPPT), 100-(horizons$CLYPPT+horizons$SLTPPT), horizons$SNDPPT)
## these are all top-soil measurements:
horizons$UHDICM <- 0
horizons$LHDICM <- 30
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2

SPROPS.OFRA <- join(horizons[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","CLYPPT","SNDPPT","SLTPPT","PHIHOX","ORCDRC","CECSUM","BLD")], SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
SPROPS.OFRA <- SPROPS.OFRA[!is.na(SPROPS.OFRA$LONWGS84) & !is.na(SPROPS.OFRA$LATWGS84) & !is.na(SPROPS.OFRA$DEPTH),]
str(SPROPS.OFRA)
summary(SPROPS.OFRA$BLD)
summary(SPROPS.OFRA$ORCDRC) ## mean = 1.7%
## 7954
save(SPROPS.OFRA, file="SPROPS.OFRA.rda")
plot(SPROPS.OFRA$LONWGS84, SPROPS.OFRA$LATWGS84, pch="+")