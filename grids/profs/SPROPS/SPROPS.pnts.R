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
## we assume - top-soil sub-soil the same
SPROPS.sim <- rbind(SPROPS.sim, SPROPS.sim)
##  http://www.jstor.org/stable/30055331
SPROPS.sim$ORCDRC = 0
SPROPS.sim$SNDPPT = 98
SPROPS.sim$SLTPPT = 2
SPROPS.sim$CLYPPT = 0
SPROPS.sim$PHIHOX = 8.1
SPROPS.sim$UHDICM = c(rep(0, nrow(deserts.pnt)), rep(5, nrow(deserts.pnt))) 
SPROPS.sim$LHDICM = c(rep(5, nrow(deserts.pnt)), rep(200, nrow(deserts.pnt)))
SPROPS.sim$DEPTH = c(rep(2.5, nrow(deserts.pnt)), rep(97.5, nrow(deserts.pnt)))
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
## Filter out remaining values outside of physical range:
all.pnts$SNDPPT <- ifelse(all.pnts$SNDPPT<0|all.pnts$SNDPPT>100, NA, all.pnts$SNDPPT)
hist(all.pnts$SNDPPT, col="blue", xlab="Sand in %", breaks=40, main=sum(!is.na(all.pnts$SNDPPT)))
all.pnts$SLTPPT <- ifelse(all.pnts$SLTPPT<0|all.pnts$SLTPPT>100, NA, all.pnts$SLTPPT)
hist(all.pnts$SLTPPT, col="blue", xlab="Silt in %", breaks=40, main=sum(!is.na(all.pnts$SLTPPT)))
all.pnts$CLYPPT <- ifelse(all.pnts$CLYPPT<0|all.pnts$CLYPPT>100, NA, all.pnts$CLYPPT)
hist(all.pnts$CLYPPT, col="blue", xlab="Clay in %", breaks=40, main=sum(!is.na(all.pnts$CLYPPT)))
all.pnts$PHIHOX <- ifelse(all.pnts$PHIHOX<2|all.pnts$PHIHOX>12, NA, all.pnts$PHIHOX)
hist(all.pnts$PHIHOX, col="blue", xlab="pH in H2O", breaks=40, main=sum(!is.na(all.pnts$PHIHOX)))
all.pnts$PHIKCL <- ifelse(all.pnts$PHIKCL<2|all.pnts$PHIKCL>12, NA, all.pnts$PHIKCL)
hist(all.pnts$PHIKCL, col="blue", xlab="pH in KCl", breaks=40, main=sum(!is.na(all.pnts$PHIKCL)))
#all.pnts$ORCDRC <- ifelse(all.pnts$ORCDRC>600, NA, all.pnts$ORCDRC)
hist(log1p(all.pnts$ORCDRC), col="blue", xlab="log-Organic carbon", breaks=40, main=sum(!is.na(all.pnts$ORCDRC)))
all.pnts$BLD <- ifelse(all.pnts$BLD<100|all.pnts$BLD>3000, NA, all.pnts$BLD)
hist(all.pnts$BLD, col="blue", xlab="Bulk density", breaks=40, main=sum(!is.na(all.pnts$BLD)))
all.pnts$CECSUM <- ifelse(all.pnts$CECSUM<0, 0, ifelse(all.pnts$CECSUM>600, NA, all.pnts$CECSUM))
hist(log1p(all.pnts$CECSUM), col="blue", xlab="log-CEC", breaks=40, main=sum(!is.na(all.pnts$CECSUM)))
summary(all.pnts)
save(all.pnts, file="all.pnts.rda")

SPROPS.pnts <- as.data.frame(all.pnts[!duplicated(all.pnts$LOC_ID),c("LOC_ID","SOURCEID","SOURCEDB","LONWGS84","LATWGS84")])
SPROPS.pnts <- SPROPS.pnts[!is.na(SPROPS.pnts$LONWGS84),]
coordinates(SPROPS.pnts) <- ~ LONWGS84+LATWGS84
proj4string(SPROPS.pnts) <- "+proj=longlat +datum=WGS84"
length(SPROPS.pnts)
summary(as.factor(SPROPS.pnts$SOURCEDB))
save(SPROPS.pnts, file="SPROPS.pnts.rda")

## world plot - overlay and plot points and maps:
s.pnts <- SPROPS.pnts[!SPROPS.pnts$SOURCEDB %in% c("NAFORMA","AU_NatSoil")&SPROPS.pnts@coords[,2]>-65&SPROPS.pnts@coords[,2]<78,]
png(file = "Fig_global_distribution_SPROPS.png", res = 150, width = 2000, height = 900)
windows(width = 20, height = 9)
dev.off()
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(s.pnts[-which(s.pnts$SOURCEDB=="Simulated"),], pch=21, bg=alpha("green", 0.6), cex=.8, col="black")
points(s.pnts[which(s.pnts$SOURCEDB=="Simulated"),], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
dev.off()