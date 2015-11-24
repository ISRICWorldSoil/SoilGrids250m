## Prepare callibration points for soil organic carbon, pH, BLD, texture fractions and similar (SoilGrids250m)
## Tom.Hengl@isric.org

library(plyr)
library(stringr)
library(sp)
library(GSIF)
library(rgdal)
library(lattice)
library(scales)
library(plotKML)
library(maps)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")

## list of input data sets:
tax.lst <- list.files(path="G:\\soilstorage\\SoilData", pattern=glob2rx("SPROPS.*.rda"), full.names=TRUE, recursive=TRUE)
tax.lst
in.lst <- lapply(tax.lst, load, .GlobalEnv)
in.lst <- lapply(in.lst, function(x){as.data.frame(get(x))})
## 24 data sets

## Simulated desert soils:
load("../TAXOUSDA/deserts.pnt.rda")
SPROPS.sim <- as.data.frame(spTransform(deserts.pnt, CRS("+proj=longlat +datum=WGS84")))
SPROPS.sim[,1] <- NULL
SPROPS.sim <- plyr::rename(SPROPS.sim, c("x"="LONWGS84", "y"="LATWGS84"))
SPROPS.sim$SOURCEID <- paste("SIM", 1:nrow(SPROPS.sim), sep="_")
SPROPS.sim$SOURCEDB = "Simulated"
##  http://www.jstor.org/stable/30055331
SPROPS.sim$ORCDRC = 0
SPROPS.sim$SNDPPT = 98
SPROPS.sim$SLTPPT = 2
SPROPS.sim$CLYPPT = 0
SPROPS.sim$PHIHOX = 8.1
SPROPS.sim$UHDICM = 0
SPROPS.sim$LHDICM = 200
SPROPS.sim$DEPTH = 100
str(SPROPS.sim)

## add to the list:
in.lst[[length(in.lst)+1]] <- SPROPS.sim
## Bind everything together:
all.pnts <- dplyr::rbind_all(in.lst)
str(all.pnts)
## 759,259
all.pnts$LOC_ID <- as.factor(paste("ID", all.pnts$LONWGS84, all.pnts$LATWGS84, sep="_"))
summary(!duplicated(all.pnts$LOC_ID))
## 146,605 unique locations!
save(all.pnts, file="all.pnts.rda")

SPROPS.pnts <- as.data.frame(all.pnts[!duplicated(all.pnts$LOC_ID),c("LOC_ID","SOURCEID","SOURCEDB","LONWGS84","LATWGS84")])
SPROPS.pnts <- SPROPS.pnts[!is.na(SPROPS.pnts$LONWGS84),]
coordinates(SPROPS.pnts) <- ~ LONWGS84+LATWGS84
proj4string(SPROPS.pnts) <- "+proj=longlat +datum=WGS84"
length(SPROPS.pnts)
summary(as.factor(SPROPS.pnts$SOURCEDB))
save(SPROPS.pnts, file="SPROPS.pnts.rda")

## world plot - overlay and plot points and maps:
no.plt <- SPROPS.pnts@coords[,2]>-65&SPROPS.pnts@coords[,2]<85
png(file = "Fig_global_distribution_SPROPS.png", res = 150, width = 2000, height = 900)
windows(width = 20, height = 9)
dev.off()
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(SPROPS.pnts[-which(SPROPS.pnts$SOURCEDB=="Simulated"&no.plt),], pch=21, bg=alpha("green", 0.6), cex=.8, col="black")
points(SPROPS.pnts[which(SPROPS.pnts$SOURCEDB=="Simulated"&no.plt),], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
dev.off()