## Soil Pedon Carbon and Nitrogen Data for Alaska: An Analysis and Update; doi:10.4236/ojss.2013.32015
## Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(plotKML)

Alaska <- read.csv("Soil_Pedon_Carbon_and_Nitrogen_Data_for_Alaska.csv")
str(Alaska)
Alaska <- rename(Alaska, c("Soil.Class"="TAXNUSDA", "Long"="LONWGS84", "Lat"="LATWGS84", "TopDepth"="UHDICM", "BtmDepth"="LHDICM", "Rock.frag...tRf."="CRFVOL", "Bulk.Density"="BLD", "SOC"="ORCDRC"))
Alaska$SOURCEID <- paste("Alaska", Alaska$General, sep="_")
Alaska$SAMPLEID <- make.unique(paste("Alaska", Alaska$General, sep="_"))

## hits on visual check:
## 341 obs with deeper 'upper' horizons than lower
View(Alaska[which(Alaska$LHDICM-Alaska$UHDICM<=0),])
#51 obs with 0 OC  
str(Alaska[which(Alaska$ORCDRC<=0),])
View(Alaska[which(Alaska$ORCDRC<=0),])

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------
 
TAXOUSDA.Alaska <- Alaska[!duplicated(Alaska$SOURCEID),c("SOURCEID","LONWGS84","LATWGS84", 
"TAXNUSDA")]
## remove missing values:
TAXOUSDA.Alaska <- TAXOUSDA.Alaska[!is.na(TAXOUSDA.Alaska$LONWGS84)&!is.na(TAXOUSDA.Alaska$TAXNUSDA)&nchar(paste(TAXOUSDA.Alaska$TAXNUSDA))>0,]
TAXOUSDA.Alaska$SOURCEDB <- "Alaska"
#View(TAXOUSDA.Alaska)
## 226 profiles!
## convert to SPDF:
coordinates(TAXOUSDA.Alaska) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.Alaska) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.Alaska["TAXNUSDA"])
save(TAXOUSDA.Alaska, file="TAXOUSDA.Alaska.rda")

##visual check 
## no detected issues 
library(tools)
View(TAXOUSDA.Alaska$TAXNUSDA)
print(showNonASCII(as.character(TAXOUSDA.Alaska$TAXNUSDA)))


# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

Alaska$ORCDRC <- as.numeric(paste(Alaska$ORCDRC)) * 10
summary(Alaska$ORCDRC) ## mean = 12%
Alaska$BLD <- as.numeric(paste(Alaska$BLD)) * 1000
summary(Alaska$BLD)
Alaska$UHDICM <- as.numeric(paste(Alaska$UHDICM))
Alaska$LHDICM <- as.numeric(paste(Alaska$LHDICM))
Alaska$SOURCEDB = "Alaska"
Alaska$DEPTH <- Alaska$UHDICM + (Alaska$LHDICM - Alaska$UHDICM)/2

SPROPS.Alaska <- Alaska[,c("SOURCEID","SAMPLEID","SOURCEDB","LONWGS84","LATWGS84","UHDICM","LHDICM","DEPTH","ORCDRC","BLD")]
SPROPS.Alaska <- SPROPS.Alaska[!is.na(SPROPS.Alaska$LONWGS84) & !is.na(SPROPS.Alaska$LATWGS84) & !is.na(SPROPS.Alaska$DEPTH),]
str(SPROPS.Alaska)
## 4062
save(SPROPS.Alaska, file="SPROPS.Alaska.rda")
plot(SPROPS.Alaska$LONWGS84, SPROPS.Alaska$LATWGS84, pch="+")
