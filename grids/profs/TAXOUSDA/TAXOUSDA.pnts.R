## Prepare callibration points for TAXOUSDA (SoilGrids250m)
## By Tom.Hengl@isric.org

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
tax.lst <- list.files(path="G:\\soilstorage\\SoilData", pattern=glob2rx("TAXOUSDA.*.rda"), full.names=TRUE, recursive=TRUE)
in.lst <- lapply(tax.lst, load, .GlobalEnv)
in.lst <- lapply(in.lst, function(x){as.data.frame(get(x))})

## Simulate desert soils:
data(landmask20km)
landmask20km <- as(landmask20km["suborder"], "SpatialPixelsDataFrame")
summary(landmask20km$suborder)
Cryids.sim <- spsample(landmask20km[landmask20km$suborder=="Cryids",], type="random", n=200)
Psamments.sim <- spsample(landmask20km[landmask20km$suborder=="Shifting Sand",], type="random", n=500)
Orthents.sim <- spsample(landmask20km[landmask20km$suborder=="Rock",], type="random", n=200)
Gelepts.sim <- spsample(landmask20km[landmask20km$suborder=="Gelepts",], type="random", n=200)
Gelands.sim <- spsample(landmask20km[landmask20km$suborder=="Gelands",], type="random", n=100)
## shoul we also simulate soils in tropics?? TH: also heavily under-represented
#Torrox.sim <- spsample(landmask20km[landmask20km$suborder=="Torrox",], type="random", n=80)
#Ustox.sim <- spsample(landmask20km[landmask20km$suborder=="Ustox",], type="random", n=80)
#Udox.sim <- spsample(landmask20km[landmask20km$suborder=="Udox",], type="random", n=80)
#Perox.sim <- spsample(landmask20km[landmask20km$suborder=="Perox",], type="random", n=80)
#Anthrepts.sim <- spsample(landmask20km[landmask20km$suborder=="Anthrepts",], type="random", n=20)
## Nicaragua points (http://hdl.handle.net/1902.1/20164)
NI.pol <- readOGR("../../../soilstorage/PolygonMaps/Nicaragua/taxonnic_wgs84.shp", "taxonnic_wgs84")
NI.sim <- spsample(NI.pol, type="random", n=150)
NI.ov <- cbind(over(NI.sim, NI.pol["SUB_GRUPO"]), NI.sim@coords)
NI.ov <- NI.ov[!is.na(NI.ov$SUB_GRUPO),]
names(NI.ov)[1] <- "TAXOUSDA"

n.lst <- c("Cryids", "Psamments", "Orthents", "Gelepts", "Gelands") ## "Torrox", "Ustox", "Udox", "Perox"
TAXOUSDA.sim <- lapply(list(Cryids.sim, Psamments.sim, Orthents.sim, Gelepts.sim, Gelands.sim), spTransform, CRS("+proj=longlat +datum=WGS84"))  ## , Torrox.sim, Ustox.sim, Udox.sim, Perox.sim
TAXOUSDA.sim.df <- list(NULL)
for(j in 1:length(TAXOUSDA.sim)){
   TAXOUSDA.sim.df[[j]] <- cbind(as.data.frame(TAXOUSDA.sim[[j]]), TAXOUSDA=rep(n.lst[j], length(TAXOUSDA.sim[[j]])))
}
TAXOUSDA.sim.df[[length(n.lst)+1]] <- NI.ov
TAXOUSDA.sim <- do.call(rbind, TAXOUSDA.sim.df)
str(TAXOUSDA.sim)
TAXOUSDA.sim <- plyr::rename(TAXOUSDA.sim, c("x"="LONWGS84", "y"="LATWGS84"))
TAXOUSDA.sim$SOURCEID <- paste("SIM", 1:nrow(TAXOUSDA.sim), sep="_")
TAXOUSDA.sim$SOURCEDB = "Simulated"
## add to the list (12):
in.lst[[length(in.lst)+1]] <- TAXOUSDA.sim

## Bind everything together:
all.pnts <- dplyr::rbind_all(in.lst)
str(all.pnts)
## 53,628
all.pnts <- as.data.frame(all.pnts)
coordinates(all.pnts) <- ~ LONWGS84+LATWGS84
proj4string(all.pnts) <- "+proj=longlat +datum=WGS84"
## Remove spatial duplicates (except for the NCSS data):
all.pnts$LOC_ID <- as.factor(paste(all.pnts@coords[,1], all.pnts@coords[,2], sep="_"))
summary(!duplicated(all.pnts$LOC_ID))
## 8813 duplicate points
selP <- all.pnts$SOURCEDB=="NCSS"
summary(selP)
#selP.pnts <- all.pnts[selP,]
sel.pnts <- all.pnts[!selP,]
## remove duplicates:
#selP.pnts <- selP.pnts[!duplicated(selP.pnts$LOC_ID),]
sel.pnts <- sel.pnts[!duplicated(sel.pnts$LOC_ID),]
TAXOUSDA.pnts <- rbind(all.pnts[selP,], sel.pnts)
length(TAXOUSDA.pnts)
## 51062
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
## FINAL NUMBER OF POINTS: 54,311 points!
summary(as.factor(TAXOUSDA.pnts$SOURCEDB))
#    CanSIS        ISIS        NCSS        WISE      AfSPDB 
#       4298         461       33501        2482        1921 
#      CIFOR      eSOTER RadamBrasil    IranSoil   Simulated 
#        212         118        6305        1535        1163 
# AU_NatSoil 
#       2076
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