## SoilGrids.org project / selection of organic soil profiles from Tropics (prapred by Nadine Herold <c9scna@gmail.com>)
## Tom.Hengl@isric.org

load(".RData")
library(plyr)
library(stringr)
library(rgdal)
library(plotKML)
library(GSIF)

## soil data out of literature:
CIFOR <- read.csv("Organic_soil_profiles_from_literature_tropics_CIFOR.csv", as.is =TRUE)
CIFOR$SOURCEID <- iconv(make.unique(gsub(" ", "_", str_trim(paste(CIFOR$author, CIFOR$year, strtrim(CIFOR$source, 8), strtrim(CIFOR$site, 8), sep="_")))), to="ASCII", sub="")
## rename:
CIFOR <- rename(CIFOR, replace=c("modelling.x"="LONWGS84", "modelling.y"="LATWGS84"))
CIFOR$SOURCEDB <- "CIFOR"
CIFOR$TIMESTRR <- as.Date(CIFOR$year, format="%Y")
summary(CIFOR$TIMESTRR)
#write.csv(CIFOR, file="CIFOR_filtered.csv")
horizons <- read.csv("CIFOR_filtered.csv")
horizons <- getHorizons(horizons, idcol="SOURCEID", sel=c("Upper", "Lower", "SOC", "BD", "Peat"))
horizons$SAMPLEID <- make.unique(as.character(horizons$SOURCEID))
horizons$UHDICM <- ifelse(is.na(horizons$Upper), 0, horizons$Upper) 
horizons$LHDICM <- ifelse(is.na(horizons$Lower), horizons$Upper+50, horizons$Lower)
horizons$ORCDRC <- horizons$SOC * 10
horizons$BLD <- horizons$BD * 1000
horizons$DEPTH <- horizons$UHDICM + (horizons$LHDICM - horizons$UHDICM)/2
## many missing values...

## hits on visual check:
## 3 LHDICM values < 800cm
str(horizons[which(horizons$LHDICM>800),])
## 4 obs with ORC higher than 600
str(horizons[which(horizons$ORCDRC>=600),])

# ------------------------------------------------------------
# export TAXONOMY DATA
# ------------------------------------------------------------

TAXOUSDA.CIFOR <- CIFOR[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TIMESTRR","TAXOUSDA")]
TAXOUSDA.CIFOR <- TAXOUSDA.CIFOR[!is.na(TAXOUSDA.CIFOR$TAXOUSDA)&nchar(TAXOUSDA.CIFOR$TAXOUSDA)>0,]
#View(TAXOUSDA.CIFOR)
## 226 profiles!
coordinates(TAXOUSDA.CIFOR) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.CIFOR) <- "+proj=longlat +datum=WGS84"
#plotKML(TAXOUSDA.CIFOR["TAXOUSDA"])
save(TAXOUSDA.CIFOR, file="TAXOUSDA.CIFOR.rda")

TAXNWRB.CIFOR <- CIFOR[,c("SOURCEID","SOURCEDB","TIMESTRR","LONWGS84","LATWGS84","TAXNWRB")]
TAXNWRB.CIFOR <- TAXNWRB.CIFOR[!is.na(TAXNWRB.CIFOR$TAXNWRB)&nchar(paste(TAXNWRB.CIFOR$TAXNWRB))>0,]
str(TAXNWRB.CIFOR)
## 247 profiles
coordinates(TAXNWRB.CIFOR) <- ~ LONWGS84+LATWGS84
proj4string(TAXNWRB.CIFOR) <- "+proj=longlat +datum=WGS84"
plotKML(TAXNWRB.CIFOR["TAXNWRB"])
save(TAXNWRB.CIFOR, file="TAXNWRB.CIFOR.rda")

##visual check - no detected issues 

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

SPROPS.CIFOR <- join(horizons[,c("SOURCEID","SAMPLEID","UHDICM","LHDICM","DEPTH","BLD","ORCDRC")], CIFOR[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TIMESTRR")], type="left")
SPROPS.CIFOR <- SPROPS.CIFOR[!is.na(SPROPS.CIFOR$LONWGS84) & !is.na(SPROPS.CIFOR$LATWGS84) & !is.na(SPROPS.CIFOR$DEPTH),]
View(SPROPS.CIFOR)
## only 840 left, but all very high values!
save(SPROPS.CIFOR, file="SPROPS.CIFOR.rda")

# ------------------------------------------------------------
# Soil organic carbon stock (kg / m2)
# ------------------------------------------------------------

SOCS.CIFOR <- CIFOR[,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB","TIMESTRR")]
SOCS.CIFOR$YEAR = CIFOR$year
SOCS.CIFOR$dSOCS_100cm = CIFOR$C.stock..MgC.ha..0.1.m/10
SOCS.CIFOR$dSOCS_200cm = CIFOR$C.stock..MgC.ha..0.2.m/10
summary(SOCS.CIFOR$dSOCS_100cm)
summary(SOCS.CIFOR$dSOCS_200cm)
coordinates(SOCS.CIFOR) <- ~ LONWGS84 + LATWGS84
proj4string(SOCS.CIFOR) = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
plot(SOCS.CIFOR)
save(SOCS.CIFOR, file="SOCS.CIFOR.rda")
