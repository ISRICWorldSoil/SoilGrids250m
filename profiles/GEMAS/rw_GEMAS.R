## Reading and writing of the GEMAS (http://gemas.geolba.ac.at/Download_GEMAS.htm);
## Reimann, C., Demetriades, A., Birke, M., Eggen, O. A., Filzmoser, P., Kriete, C. & EuroGeoSurveys Geochemistry Expert Group, 2012. The EuroGeoSurveys Geochemical Mapping of Agricultural and grazing land Soils project (GEMAS) - Evaluation of quality control results of particle size estimation by MIR® prediction, Pb-isotope and MMI® extraction analyses and results of the GEMAS ring test for the standards Ap and Gr. NGU Report 2012.051, 136 pp.

setwd("/data/soilstorage/SoilData/GEMAS")
load(".RData")
library(aqp)
library(GSIF)
library(maptools)
#library(gdata)
library(plyr)
profs <- read.csv("GEMAS.csv", sep=",", header = TRUE)
str(profs)

## ID column:
plyr:::nunique(profs$ID)

## rename columns:
profs <- plyr::rename(profs, c("clay"="CLYPPT", "silt"="SLTPPT", "pH_CaCl2"="PHICAL", "TOC"="ORCDRC", "YCOO"="LATWGS84", "XCOO"="LONWGS84", "CEC"="CECSUM"))
profs$SOURCEID = paste("GEMAS",profs$ID,sep="_")
str(profs)

## check / convert values where necessary
## Organic carbon in permilles:
profs$ORCDRC = profs$ORCDRC*10
summary(profs$ORCDRC)
summary(profs$PHICAL)
profs$SNDPPT = 100 - (profs$SLTPPT + profs$CLYPPT)
hist(profs$SNDPPT)

profs$TIMESTRR <- as.Date("2008", format="%Y")
## subset to complete data:
sel.c <- !is.na(profs$LONWGS84)&!is.na(profs$LATWGS84)
profs.f <- profs[sel.c,]

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

profs.f$SAMPLEID <- make.unique(paste0("GEMAS", profs.f$SOURCEID))
profs.f$SOURCEDB = "GEMAS"
profs.f$DEPTH <- profs.f$UHDICM + (profs.f$LHDICM - profs.f$UHDICM)/2
SPROPS.GEMAS <- profs.f[,c("SOURCEID","SAMPLEID","SOURCEDB","LONWGS84","LATWGS84","TIMESTRR","UHDICM","LHDICM","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHICAL","ORCDRC","CECSUM")]
str(SPROPS.GEMAS)
## 4131
save(SPROPS.GEMAS, file="SPROPS.GEMAS.rda")
plot(SPROPS.GEMAS$LONWGS84, SPROPS.GEMAS$LATWGS84, pch="+")
saveRDS(profs.f[,c("SOURCEID","SAMPLEID","SOURCEDB","LONWGS84","LATWGS84","TIMESTRR","UHDICM","LHDICM","DEPTH","SNDPPT","CLYPPT","SLTPPT","PHICAL","ORCDRC","CECSUM","As", "Cd", "Cu", "Pb", "Zn")], file="pnts_GEMAS.rds")
save.image()

