## A new data set for estimating organic carbon storage to 3m depth in soils of the northern circumpolar permafrost region; doi:10.5194/essd-5-393-2013
## Download from: http://bolin.su.se/data/ncscd/pedon.php
## Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(plotKML)

artic <- read.csv("Harden_etal_2012_Hugelius_etal_2013_cleaned_data_for_ISRIC.csv", stringsAsFactors=FALSE)
str(artic)
artic <- rename(artic, c("Suborder"="TAXOUSDA", "Long"="LONWGS84", "Lat"="LATWGS84", "Basal.Depth.cm"="UHDICM", "bulk.density.g.cm.3"="BLD", "X.C"="ORCDRC"))
artic$LHDICM <- artic$UHDICM + artic$Layer.thickness.cm
artic$SOURCEID <- paste("ART", artic$Profile.ID, sep="_")
artic$SAMPLEID <- make.unique(paste("ART", artic$Profile.ID, artic$Horizon.type, sep="_"))

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------
 
TAXOUSDA.artic <- artic[!duplicated(artic$SOURCEID),c("SOURCEID","LONWGS84","LATWGS84","TAXOUSDA")]
## rename:
TAXOUSDA.artic <- TAXOUSDA.artic[!is.na(TAXOUSDA.artic$LONWGS84)&!is.na(TAXOUSDA.artic$TAXOUSDA)&nchar(paste(TAXOUSDA.artic$TAXOUSDA))>0,]
TAXOUSDA.artic$SOURCEDB <- "Artic"
#View(TAXOUSDA.artic)
## 1175 profiles!
## convert to SPDF:
coordinates(TAXOUSDA.artic) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.artic) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.artic["TAXOUSDA"])
save(TAXOUSDA.artic, file="TAXOUSDA.artic.rda")

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

artic$ORCDRC <- as.numeric(paste(artic$ORCDRC)) * 10
summary(artic$ORCDRC)
artic$BLD <- as.numeric(paste(artic$BLD)) * 1000
summary(artic$BLD)
artic$BLD <- ifelse(artic$BLD < 100 | artic$BLD > 3000, NA, artic$BLD)
artic$SOURCEDB = "Artic"
artic$DEPTH <- artic$UHDICM + (artic$LHDICM - artic$UHDICM)/2
artic$SOURCEDB <- "Artic"

SPROPS.artic <- artic[,c("SOURCEID","SAMPLEID","SOURCEDB","LONWGS84","LATWGS84","UHDICM","LHDICM","DEPTH","ORCDRC","BLD")]
SPROPS.artic <- SPROPS.artic[!is.na(SPROPS.artic$LONWGS84) & !is.na(SPROPS.artic$LATWGS84) & !is.na(SPROPS.artic$DEPTH),]
str(SPROPS.artic)
## 7457
## Needs to be checked! the histograms show some 5000 measurements >40% of SOC
## Most likely the authors have over-represented histosols / peatlands
save(SPROPS.artic, file="SPROPS.artic.rda")
plot(SPROPS.artic$LONWGS84, SPROPS.artic$LATWGS84, pch="+")
