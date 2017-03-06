## Mangrove soil profiles from literature (prepared by TNC)  
## This data has not yet been published and license has not yet been set!

library(aqp)
library(GSIF)
library(sp)
library(maptools)
library(plyr)

profs <- read.csv("mangrove_soc_database_v10_sites.csv")
#str(profs)
profs.f <- rename(profs, c("Site.name"="SOURCEID", "Longitude"="LONWGS84", "Latitude"="LATWGS84"))
profs.f$TIMESTRR <- as.Date(profs.f$Years_collected, format="%Y")
profs.f$SOURCEDB = "MangrovesDB"
profs.f$LONWGS84 = as.numeric(paste(profs.f$LONWGS84))
profs.f$LATWGS84 = as.numeric(paste(profs.f$LATWGS84))

hors <- read.csv("mangrove_soc_database_v10_horizons.csv")
#str(hors)
hors.f <- rename(hors, c("Site.name"="SOURCEID", "U_depth"="UHDICM", "L_depth"="LHDICM", "OC_final"="ORCDRC", "BD_final"="BLD"))
summary(hors.f$ORCDRC)
hors.f$ORCDRC <- hors.f$ORCDRC*10
summary(hors.f$BLD)
hors.f$BLD <- hors.f$BLD*1000
## convert depths to cm:
summary(hors.f$UHDICM)
hors.f$UHDICM = hors.f$UHDICM*100
summary(hors.f$LHDICM)
hors.f$LHDICM = hors.f$LHDICM*100

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

hors.f$DEPTH <- hors.f$UHDICM + (hors.f$LHDICM - hors.f$UHDICM)/2
hist(hors.f$DEPTH)
SPROPS.MangrovesDB <- join(profs.f[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84")], hors.f[,c("SOURCEID","UHDICM","LHDICM","DEPTH","BLD","ORCDRC")])
SPROPS.MangrovesDB = SPROPS.MangrovesDB[!is.na(SPROPS.MangrovesDB$LONWGS84) & !is.na(SPROPS.MangrovesDB$LATWGS84) & !is.na(SPROPS.MangrovesDB$DEPTH),]
str(SPROPS.MangrovesDB)
## 7901
save(SPROPS.MangrovesDB, file="SPROPS.MangrovesDB.rda")
save.image()
plot(SPROPS.MangrovesDB$LONWGS84, SPROPS.MangrovesDB$LATWGS84, pch="+")
