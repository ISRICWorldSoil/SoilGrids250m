## Prepare callibration points for mapping soil depth (SoilGrids250m)
## Tom.Hengl@isric.org / Wei.Shangguan
## 250 indicates =>250

library(plyr)
library(stringr)
library(sp)
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
tax.lst <- list.files(path="G:\\soilstorage\\SoilData", pattern=glob2rx("BDR.*.rda"), full.names=TRUE, recursive=TRUE)
tax.lst
in.lst <- lapply(tax.lst, load, .GlobalEnv)
in.lst <- lapply(in.lst, function(x){as.data.frame(get(x))})
## 15 data sets

## Well data:
load("wells.depth.rda")
str(wells.depth) ## 1.5M points
BDR.wells <- as.data.frame(wells.depth)
BDR.wells$BDTICM <- BDR.wells$BDRICM ## absolute depth to R
hist(log1p(BDR.wells$BDRICM)) ## log-normal distribution
summary(BDR.wells$BDRICM>250) ## 380k points shallow
BDR.wells$SOURCEDB = "Wells"
nrow(BDR.wells)
## 1.5M points

## Simulated desert soils:
load("../TAXOUSDA/deserts.pnt.rda")
load("../TAXOUSDA/barerock.pnt.rda")
d <- as.data.frame(spTransform(deserts.pnt, CRS("+proj=longlat +datum=WGS84")))
d <- plyr::rename(d, c("desertPR_sin"="BDRICM", "x"="LONWGS84", "y"="LATWGS84"))
d$BDRICM <- 250
b <- as.data.frame(spTransform(barerock.pnt, CRS("+proj=longlat +datum=WGS84")))
b <- plyr::rename(b, c("barerockPR_sin"="BDRICM", "x"="LONWGS84", "y"="LATWGS84"))
b$BDRICM <- 5
BDR.sim <- rbind(d,b)
BDR.sim$SOURCEID <- paste("SIM", 1:nrow(BDR.sim), sep="_")
BDR.sim$SOURCEDB = "Simulated"
str(BDR.sim) ## 541 point

## add to the list:
in.lst[[length(in.lst)+1]] <- BDR.sim
in.lst[[length(in.lst)+1]] <- BDR.wells
## Bind everything together:
BDR_all.pnts <- dplyr::rbind_all(in.lst)
str(BDR_all.pnts)
## 1,7M points
BDR_all.pnts$LOC_ID <- as.factor(paste("ID", BDR_all.pnts$LONWGS84, BDR_all.pnts$LATWGS84, sep="_"))
summary(!duplicated(BDR_all.pnts$LOC_ID))
## 1,5M unique locations!
summary(BDR_all.pnts$BDRICM)
## copy values from soil profile data:
BDR_all.pnts$BDTICM <- ifelse(is.na(BDR_all.pnts$BDTICM), ifelse(!BDR_all.pnts$BDRICM==250, BDR_all.pnts$BDRICM, NA),  BDR_all.pnts$BDTICM)
hist(log1p(BDR_all.pnts$BDTICM))
BDR_all.pnts$BDRICM <- ifelse(BDR_all.pnts$BDRICM>=250, 250, BDR_all.pnts$BDRICM)
BDR_all.pnts$BDRLOG <- ifelse(BDR_all.pnts$BDRICM<250, 1, 0)
str(BDR_all.pnts)
summary(as.factor(BDR_all.pnts$BDRLOG))
#      0       1 
#1,319,148  387,161
save(BDR_all.pnts, file="BDR_all.pnts.rda")
BDR_all.pnts <- as.data.frame(BDR_all.pnts)
BDR_all.pnts[1,]
BDR_all.pnts[800000,]
## BDRLOG = 0
## BDTICM = 609
## BDRICM = 250

BDR.pnts <- BDR_all.pnts[!duplicated(BDR_all.pnts$LOC_ID),c("LOC_ID","SOURCEID","SOURCEDB","LONWGS84","LATWGS84")]
BDR.pnts <- BDR.pnts[!is.na(BDR.pnts$LONWGS84),]
coordinates(BDR.pnts) <- ~ LONWGS84+LATWGS84
proj4string(BDR.pnts) <- "+proj=longlat +datum=WGS84"
length(BDR.pnts) ## 1.5M points
summary(as.factor(BDR.pnts$SOURCEDB))
save(BDR.pnts, file="BDR.pnts.rda")

## world plot - overlay and plot points and maps:
no.plt <- BDR.pnts@coords[,2]>-65&BDR.pnts@coords[,2]<85
png(file = "Fig_global_distribution_BDR.png", res = 150, width = 2000, height = 900)
windows(width = 20, height = 9)
dev.off()
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
## profile data:
points(BDR.pnts[-which(BDR.pnts$SOURCEDB=="Wells"|BDR.pnts$SOURCEDB=="Simulated"&!no.plt),], pch=21, bg=alpha("red", 0.6), cex=.8, col="black")
## Wells data
points(BDR.pnts[-which(!BDR.pnts$SOURCEDB=="Wells"|BDR.pnts$SOURCEDB=="Simulated"&no.plt),], pch=21, bg=alpha("blue", 0.6), cex=.8, col="black")
points(BDR.pnts[which(BDR.pnts$SOURCEDB=="Simulated"&no.plt),], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
dev.off()