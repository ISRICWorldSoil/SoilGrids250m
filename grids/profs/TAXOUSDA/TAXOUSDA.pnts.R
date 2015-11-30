## Prepare callibration points for TAXOUSDA (SoilGrids250m)
## By Tom.Hengl@isric.org

library(plyr)
library(stringr)
library(sp)
library(GSIF)
library(rgdal)
library(raster)
library(lattice)
library(scales)
library(plotKML)
library(maps)
country.m <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country.m$names, ":"), function(x) x[1])
require(maptools)
country <- as(map2SpatialPolygons(country.m, IDs=IDs), "SpatialLines")

## list of input data sets:
tax.lst <- list.files(path="G:\\soilstorage\\SoilData", pattern=glob2rx("TAXOUSDA.*.rda"), full.names=TRUE, recursive=TRUE)
tax.lst ## 14
in.lst <- lapply(tax.lst, load, .GlobalEnv)
in.lst <- lapply(in.lst, function(x){as.data.frame(get(x))})

## Desert / barerock soils (see "Encyclopedia of Soil Science" Vol.1 P520)
deserts <- raster("../../masks/desertPR_sin.tif")
barerock <- raster("../../masks/barerockPR_sin.tif")
#deserts.pnt <- sampleRandom(deserts, size=2000, sp=TRUE)
#plot(deserts.pnt) ## ca 400 points
#save(deserts.pnt, file="deserts.pnt.rda")
#barerock.pnt <- sampleRandom(barerock, size=30000, sp=TRUE)
#plot(barerock.pnt) ## ca 200 points 
#save(barerock.pnt, file="barerock.pnt.rda")
load("deserts.pnt.rda")
load("barerock.pnt.rda")

## Nicaragua points (http://hdl.handle.net/1902.1/20164)
NI.pol <- readOGR("../../../soilstorage/PolygonMaps/Nicaragua/taxonnic_wgs84.shp", "taxonnic_wgs84")
NI.sim <- spsample(NI.pol, type="random", n=150)
NI.ov <- cbind(over(NI.sim, NI.pol["SUB_GRUPO"]), NI.sim@coords)
NI.ov <- NI.ov[!is.na(NI.ov$SUB_GRUPO),]
names(NI.ov)[1] <- "TAXOUSDA"

n.lst <- c("Psamments", "Orthents")
TAXOUSDA.sim <- lapply(list(deserts.pnt, barerock.pnt), spTransform, CRS("+proj=longlat +datum=WGS84")) 
TAXOUSDA.sim.df <- list(NULL)
for(j in 1:length(TAXOUSDA.sim)){
   TAXOUSDA.sim.df[[j]] <- cbind(as.data.frame(TAXOUSDA.sim[[j]]), TAXOUSDA=rep(n.lst[j], length(TAXOUSDA.sim[[j]])))
   TAXOUSDA.sim.df[[j]][,1] <- NULL
}
TAXOUSDA.sim.df[[length(n.lst)+1]] <- NI.ov
TAXOUSDA.sim <- do.call(rbind, TAXOUSDA.sim.df)
str(TAXOUSDA.sim)
TAXOUSDA.sim <- plyr::rename(TAXOUSDA.sim, c("x"="LONWGS84", "y"="LATWGS84"))
TAXOUSDA.sim$SOURCEID <- paste("SIM", 1:nrow(TAXOUSDA.sim), sep="_")
TAXOUSDA.sim$SOURCEDB = "Simulated"
## add simulated points to the list:
in.lst[[length(in.lst)+1]] <- TAXOUSDA.sim

## Bind everything together:
all.pnts <- dplyr::rbind_all(in.lst)
str(all.pnts)
## 55,196
all.pnts <- as.data.frame(all.pnts)
coordinates(all.pnts) <- ~ LONWGS84+LATWGS84
proj4string(all.pnts) <- "+proj=longlat +datum=WGS84"
## Remove spatial duplicates (except for the NCSS data):
all.pnts$LOC_ID <- as.factor(paste(all.pnts@coords[,1], all.pnts@coords[,2], sep="_"))
summary(!duplicated(all.pnts$LOC_ID))
## 10,347 duplicate points!
selP <- all.pnts$SOURCEDB=="NCSS"
summary(selP)
#selP.pnts <- all.pnts[selP,]
sel.pnts <- all.pnts[!selP,]
## remove duplicates but keep all NCSS points:
#selP.pnts <- selP.pnts[!duplicated(selP.pnts$LOC_ID),]
sel.pnts <- sel.pnts[!duplicated(sel.pnts$LOC_ID),]
TAXOUSDA.pnts <- rbind(all.pnts[selP,], sel.pnts)
length(TAXOUSDA.pnts)
## 51580
levels(as.factor(all.pnts$SOURCEDB))

## clean-up names:
USDA_levs <- read.csv("USDA_levs.csv", na.strings = c("NA",""))
str(USDA_levs)
levs <- levels(USDA_levs$Group)[-which(levels(USDA_levs$Group) %in% c("Shifting Sand", "Rock", "Ocean", "Ice"))]
# 70 suborders
## One by one suborder name --> try to located them in the raw names
tax.lst <- list(NULL)
for(j in 1:length(levs)){
  ## remove "s" if at the end of the class name:
  pat <- ifelse(substr(levs[j], nchar(levs[j]), nchar(levs[j]))=="s", substr(levs[j], 1, nchar(levs[j])-1), levs[j])
  sel1 <- grep(pat, paste(TAXOUSDA.pnts$TAXNUSDA), ignore.case=TRUE)
  sel2 <- grep(pat, paste(TAXOUSDA.pnts$TAXOUSDA), ignore.case=TRUE)
  sel <- unique(c(sel1, sel2))
  ## there can be multiple soil orders at the same location!
  if(length(sel)>0){
    tax.lst[[j]] <- data.frame(TAXOUSDA.pnts[sel,])
    tax.lst[[j]]$TAXOUSDA.f = levs[j]
  }
}
TAXOUSDA.pnts <- do.call(rbind, tax.lst)
TAXOUSDA.pnts$TAXOUSDA.f <- as.factor(TAXOUSDA.pnts$TAXOUSDA.f)
summary(TAXOUSDA.pnts$TAXOUSDA.f)
length(TAXOUSDA.pnts$TAXOUSDA.f)
## FINAL NUMBER OF POINTS: 55,932 points!
summary(as.factor(TAXOUSDA.pnts$SOURCEDB))
#          CanSIS                                      ISIS 
#             4298                                       461 
#             NCSS                              Russia_EGRPR 
#            33501                                       623 
#       USGS_Buell                                      WISE 
#             1060                                      2481 
#           AfSPDB                                     CIFOR 
#             1921                                       212 
#           eSOTER                               RadamBrasil 
#              118                                      6305 
#           Alaska                                  IranSoil 
#              387                                      1535 
#        Simulated                                AU_NatSoil 
#              673                                      2076 
#            Artic University_of_Michigan_Biological_Station 
#              275                                         1
coordinates(TAXOUSDA.pnts) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.pnts) <- "+proj=longlat +datum=WGS84"
plotKML(TAXOUSDA.pnts["TAXOUSDA.f"], file.name="TAXOUSDA_Oct_27_2015.kml", kmz=TRUE)
str(TAXOUSDA.pnts@data)
TAXOUSDA.pnts <- TAXOUSDA.pnts[c("SOURCEDB","SOURCEID","TAXOUSDA.f")]
save(TAXOUSDA.pnts, file="TAXOUSDA.pnts.rda")

## world plot - overlay and plot points and maps:
no.plt <- TAXOUSDA.pnts@coords[,2]>-65&TAXOUSDA.pnts@coords[,2]<85
png(file = "Fig_global_distribution_TAXOUSDA.png", res = 150, width = 2000, height = 900)
windows(width = 20, height = 9)
dev.off()
par(mar=c(0,0,0,0), oma=c(0,0,0,0))
plot(country, col="darkgrey", ylim=c(-60, 85))
points(TAXOUSDA.pnts[!TAXOUSDA.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("blue", 0.6), cex=.8, col="black")
points(TAXOUSDA.pnts[TAXOUSDA.pnts$SOURCEDB=="Simulated"&no.plt,], pch=21, bg=alpha("yellow", 0.6), cex=.6, col="black")
dev.off()