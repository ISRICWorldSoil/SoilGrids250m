## Data from (http://www.iransoil.com/webgis.html)
## Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(plotKML)

Iransoil <- read.csv("Iransoil.csv")
str(Iransoil)
Iransoil <- rename(Iransoil, c("classification.soil.taxonomy."="TAXNUSDA", "longitude"="LONWGS84", "latitude"="LATWGS84", "X..sand."="SNDPPT", "X..silt."="SLTPPT", "X..clay."="CLYPPT", "pH"="PHIHOX", "Bulk.Density..b..g.cm3."="BLD", "X.O.C"="ORCDRC"))
Iransoil$SOURCEID <- make.unique(paste("IranSoil", Iransoil$LONWGS84, Iransoil$LATWGS84, sep="_"))
TAXOUSDA.Iransoil <- Iransoil[,c("SOURCEID","LONWGS84","LATWGS84", "TAXNUSDA")]
## rename:
TAXOUSDA.Iransoil <- TAXOUSDA.Iransoil[!is.na(TAXOUSDA.Iransoil$LONWGS84)&!is.na(TAXOUSDA.Iransoil$TAXNUSDA)&nchar(paste(TAXOUSDA.Iransoil$TAXNUSDA))>0,]
TAXOUSDA.Iransoil$SOURCEDB <- "IranSoil"
#View(TAXOUSDA.Iransoil)
## 226 profiles!
## convert to SPDF:
coordinates(TAXOUSDA.Iransoil) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.Iransoil) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.Iransoil["TAXNUSDA"])
save(TAXOUSDA.Iransoil, file="TAXOUSDA.Iransoil.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

## top soil only
Iransoil$UHDICM = 0
Iransoil$LHDICM = 30
Iransoil$ORCDRC <- as.numeric(paste(Iransoil$ORCDRC)) * 10
summary(Iransoil$ORCDRC)
Iransoil$BLD <- Iransoil$BLD * 1000
summary(Iransoil$BLD)
Iransoil$PHIHOX <- as.numeric(paste(Iransoil$PHIHOX))
summary(Iransoil$PHIHOX)
Iransoil$SAMPLEID <- make.unique(paste(Iransoil$SOURCEID, 1, sep="_"))
Iransoil$SOURCEDB = "IranSoil"
Iransoil$DEPTH <- Iransoil$UHDICM + (Iransoil$LHDICM - Iransoil$UHDICM)/2

SPROPS.Iransoil <- Iransoil[,c("SOURCEID","SAMPLEID","SOURCEDB","LONWGS84","LATWGS84","UHDICM","LHDICM","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHIHOX","ORCDRC","BLD")]
SPROPS.Iransoil <- SPROPS.Iransoil[!is.na(SPROPS.Iransoil$LONWGS84) & !is.na(SPROPS.Iransoil$LATWGS84) & !is.na(SPROPS.Iransoil$DEPTH),]
str(SPROPS.Iransoil)
## 2026
save(SPROPS.Iransoil, file="SPROPS.Iransoil.rda")
plot(SPROPS.Iransoil$LONWGS84, SPROPS.Iransoil$LATWGS84, pch="+")