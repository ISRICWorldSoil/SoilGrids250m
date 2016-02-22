# reference     : bodemdatabank Aardewerk-Vlaanderen-2010 [http://www.sadl.kuleuven.be/]
# producer      : by Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(sp)
library(plotKML)

SITE <- read.csv("Aardewerk-Vlaanderen-2010_Profiel.csv")
str(SITE)
## 7954 layers
s <- summary(SITE$Bodemgroep, maxsum = 140)
soiltype <- data.frame(soil_name=attr(s, "names"), count=s)
write.csv(soiltype, "Bodemgroep_count.csv")
SITE$SOURCEDB = "Aardewerk-Vlaanderen-2010"
SITE$SOURCEID <- paste(SITE$ID)
## Coordinate system: EPSG:31300
SITE.s <- SITE[!duplicated(SITE$SOURCEID)& !is.na(SITE$Coordinaat_Lambert72_X),c("SOURCEID","SOURCEDB","Coordinaat_Lambert72_X","Coordinaat_Lambert72_Y","Bodemgroep","Profilering_Datum")]
coordinates(SITE.s) <- ~ Coordinaat_Lambert72_X + Coordinaat_Lambert72_Y
proj4string(SITE.s) <- "+init=epsg:31300"
#plotKML(SITE.s[1])
#legFAO_90 <- read.csv("cleanup_Vlaanderen.csv")
SITE.s <- as.data.frame(spTransform(SITE.s, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")))
#SITE.s$TAXNWRB <- join(SITE.s, legFAO_90, type="left")$Name
SITE.s$TIMESTRR <- as.Date(paste(SITE.s$Profilering_Datum), format="%d-%m-%Y")
SITE.s <- rename(SITE.s, c("Coordinaat_Lambert72_X"="LONWGS84", "Coordinaat_Lambert72_Y"="LATWGS84"))
View(SITE.s)

horizons <- read.csv("Aardewerk-Vlaanderen-2010_Horizont.csv")
str(horizons)
horizons$ORCDRC <- horizons$Humus*10 /1.724     ## OC in permilles
hist(horizons$ORCDRC)
## few typoes:
horizons$ORCDRC <- ifelse(horizons$ORCDRC>600, NA, horizons$ORCDRC)
horizons$SAMPLEID <- make.unique(paste(horizons$Profiel_ID, horizons$Hor_nr, sep="_"))
horizons <- rename(horizons, c("Profiel_ID"="SOURCEID", "pH_KCl"="PHIKCL", "pH_H2O"="PHIHOX", "Tgroter_dan_2000"="CRFVOL", "T0_2"="CLYPPT", "Sorptiecapaciteit_Totaal"="CECSUM"))
horizons$UHDICM <- rowSums(horizons[,c("Diepte_grens_boven1", "Diepte_grens_boven2")], na.rm=TRUE)/2
horizons$LHDICM <- rowSums(horizons[,c("Diepte_grens_onder1","Diepte_grens_onder2")], na.rm=TRUE)/2
summary(horizons$UHDICM)
horizons$SNDPPT <- horizons$T50_100 + horizons$T100_200 + horizons$T200_500 + horizons$T500_1000 + horizons$T1000_2000
horizons$SLTPPT <- horizons$T2_10 + horizons$T10_20 + horizons$T20_50
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
horizons$SLTPPT <- ifelse(is.na(horizons$SLTPPT), 100-(horizons$CLYPPT+horizons$SNDPPT), horizons$SLTPPT)
horizons$SNDPPT <- ifelse(is.na(horizons$SNDPPT), 100-(horizons$CLYPPT+horizons$SLTPPT), horizons$SNDPPT)
sumTex <- rowSums(horizons[,c("SLTPPT","CLYPPT","SNDPPT")])
hist(sumTex)

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXNWRB.Vlaanderen <- SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.Vlaanderen <- TAXNWRB.Vlaanderen[!is.na(TAXNWRB.Vlaanderen$TAXNWRB)&!is.na(TAXNWRB.Vlaanderen$LONWGS84)&nchar(paste(TAXNWRB.Vlaanderen$TAXNWRB))>0,]
str(TAXNWRB.Vlaanderen)
## 101 profiles
coordinates(TAXNWRB.Vlaanderen) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.Vlaanderen) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.Vlaanderen["TAXNWRB"])
save(TAXNWRB.Vlaanderen, file="TAXNWRB.Vlaanderen.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.Vlaanderen <- join(horizons[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIHOX","PHIKCL","ORCDRC","CECSUM","CRFVOL")], SITE.s[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84")], type="left")
SPROPS.Vlaanderen <- SPROPS.Vlaanderen[!is.na(SPROPS.Vlaanderen$LONWGS84) & !is.na(SPROPS.Vlaanderen$LATWGS84) & !is.na(SPROPS.Vlaanderen$DEPTH),]
str(SPROPS.Vlaanderen)
## 41,789
save(SPROPS.Vlaanderen, file="SPROPS.Vlaanderen.rda")
plot(SPROPS.Vlaanderen$LONWGS84, SPROPS.Vlaanderen$LATWGS84, pch="+")

# ------------------------------------------------------------
# Depth to bedrock
# ------------------------------------------------------------

