## Reading and writing of the LUCAS soil profile data;
## The LUCAS_SOIL data (19,969) are the property of the European Union, represented by the European Commission, represented by the Directorate General-Joint Research Centre [http://eusoils.jrc.ec.europa.eu/projects/Lucas/Data.html];
## T?th, G., Jones, A., Montanarella, L. (eds.) 2013. LUCAS Topsoil Survey. Methodology, data and results. JRC Technical Reports. Luxembourg. Publications Office of the European Union, EUR 26102 ? Scientific and Technical Research series ? ISSN 1831-9424 (online); ISBN 978-92-79-32542-7; doi: 10.2788/97922

load(".RData")
library(aqp)
library(GSIF)
library(maptools)
library(gdata)
library(plyr)
perl <- gdata:::findPerl("perl")
installXLSXsupport(perl=perl)
tmp <- read.xls("LUCAS_TOPSOIL_v1.xlsx", perl=perl, sheet=1)
## takes 2-3 mins to read!
profs <- tmp

## ID column:
profs$SOURCEID <- as.factor(ifelse(profs$POINT_ID=="NE", paste("X", row(profs), sep=""), paste("ID", profs$POINT_ID, sep="")))
# check if there are duplicates:
length(levels(profs$SOURCEID))
plyr:::nunique(profs$SOURCEID)

## rename columns:
profs <- plyr::rename(profs, c("coarse"="CRFVOL", "clay"="CLYPPT", "silt"="SLTPPT", "sand"="SNDPPT", "pH_in_H2O"="PHIHOX", "pH_in_CaCl"="PHICAL", "OC"="ORCDRC", "GPS_LAT"="LATWGS84",	"GPS_LONG"="LONWGS84", "CEC"="CECSUM"))
str(profs)

## check / convert values where necessary
## Organic carbon in permilles:
summary(profs$ORCDRC)
summary(profs$PHIHOX)
summary(profs$PHIKCL)
summary(rowSums(profs[,c("CLYPPT", "SLTPPT", "SNDPPT")]) > 102 | rowSums(profs[,c("CLYPPT", "SLTPPT", "SNDPPT")]) < 98)
## 1128 significantly different from 100%
profs$UHDICM <- 0
profs$LHDICM <- 20
profs$TIMESTRR <- as.Date("2009", format="%Y")
## subset to complete data:
sel.c <- !is.na(profs$LONWGS84)&!is.na(profs$LATWGS84)
profs.f <- profs[sel.c,]

# ------------------------------------------------------------
# All soil properties
# ------------------------------------------------------------

profs.f$SAMPLEID <- make.unique(paste(profs.f$sample_ID))
profs.f$SOURCEDB = "LUCAS"
profs.f$DEPTH <- profs.f$UHDICM + (profs.f$LHDICM - profs.f$UHDICM)/2
SPROPS.LUCAS <- profs.f[,c("SOURCEID","SAMPLEID","SOURCEDB","LONWGS84","LATWGS84","TIMESTRR","UHDICM","LHDICM","DEPTH","SNDPPT","CRFVOL","CLYPPT","SLTPPT","PHIHOX","PHICAL","ORCDRC","CECSUM")]
str(SPROPS.LUCAS)
## 19,899
save(SPROPS.LUCAS, file="SPROPS.LUCAS.rda")
plot(SPROPS.LUCAS$LONWGS84, SPROPS.LUCAS$LATWGS84, pch="+")