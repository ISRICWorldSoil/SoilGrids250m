## World compilation of USDA suborders taxa
## by: Tom.Hengl@isric.org

library(aqp)
library(plyr)
library(stringr)
library(sp)
library(rgdal)
library(plotKML)
library(GSIF)
library(dplyr)
library(snowfall)
library(rgdal)

tax.lst <- c("USDA/NCSS/TAXOUSDA.NCSS.rda", "IranSoil/TAXOUSDA.Iransoil.rda", "ISCN/TAXOUSDA.ISCN.rda", "CIFOR/TAXOUSDA.CIFOR.rda", "WISE/TAXOUSDA.WISE.rda", "ISIS/TAXOUSDA.ISIS.rda")
in.lst <- lapply(paste0("../../SoilData/", tax.lst), load, .GlobalEnv)
in.lst <- lapply(in.lst, function(x){as.data.frame(get(x))})
col.legend <- read.csv("TAXOUSDA_legend.csv")
col.legend <- col.legend[!is.na(col.legend$R),]
col.legend$COLOR <- rgb(red=col.legend$R/255, green=col.legend$G/255, blue=col.legend$B/255)

## desert soils:
data(landmask20km)
landmask20km <- as(landmask20km["suborder"], "SpatialPixelsDataFrame")
summary(landmask20km$suborder)
Cryids.sim <- spsample(landmask20km[landmask20km$suborder=="Cryids",], type="random", n=200)
Psamments.sim <- spsample(landmask20km[landmask20km$suborder=="Shifting Sand",], type="random", n=500)
Orthents.sim <- spsample(landmask20km[landmask20km$suborder=="Rock",], type="random", n=200)
Gelepts.sim <- spsample(landmask20km[landmask20km$suborder=="Gelepts",], type="random", n=200)
#Anthrepts.sim <- spsample(landmask20km[landmask20km$suborder=="Anthrepts",], type="random", n=20)

n.lst <- c("Cryids", "Psamments", "Orthents", "Gelepts")
TAXOUSDA.sim <- lapply(list(Cryids.sim, Psamments.sim, Orthents.sim, Gelepts.sim), spTransform, CRS("+proj=longlat +datum=WGS84"))
TAXOUSDA.sim.df <- list(NULL)
for(j in 1:length(TAXOUSDA.sim)){
   TAXOUSDA.sim.df[[j]] <- cbind(as.data.frame(TAXOUSDA.sim[[j]]), TAXOUSDA=rep(n.lst[j], length(TAXOUSDA.sim[[j]])))
}
TAXOUSDA.sim <- do.call(rbind, TAXOUSDA.sim.df)
str(TAXOUSDA.sim)
TAXOUSDA.sim <- plyr::rename(TAXOUSDA.sim, c("x"="LONWGS84", "y"="LATWGS84"))
TAXOUSDA.sim$SOURCEID <- paste("SIM", 1:nrow(TAXOUSDA.sim), sep="_")
TAXOUSDA.sim$SOURCEDB = "Simulated"
in.lst[[7]] <- TAXOUSDA.sim

## takes few secs:
TAXOUSDA.pnts <- dplyr::rbind_all(in.lst)
str(TAXOUSDA.pnts)
TAXOUSDA.pnts <- as.data.frame(TAXOUSDA.pnts)
coordinates(TAXOUSDA.pnts) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.pnts) <- "+proj=longlat +datum=WGS84"
## clean-up names:
USDA_levs <- read.csv("USDA_levs.csv", na.strings = c("NA",""))
str(USDA_levs)
levs <- levels(USDA_levs$Group)[-which(levels(USDA_levs$Group) %in% c("Shifting Sand", "Rock", "Ocean", "Ice"))]
TAXOUSDA.pnts$TAXOUSDA.f <- NA
## Go one by one suborder name and try to located them in the raw names
for(j in 1:length(levs)){
  pat <- ifelse(substr(levs[j], nchar(levs[j]), nchar(levs[j]))=="s", substr(levs[j], 1, nchar(levs[j])-1), nchar(levs[j]))
  sel1 <- grep(pat, paste(TAXOUSDA.pnts$TAXNUSDA), ignore.case=TRUE)
  sel2 <- grep(pat, paste(TAXOUSDA.pnts$TAXOUSDA), ignore.case=TRUE)
  sel <- unique(c(sel1, sel2))
  if(length(sel)>0){
    TAXOUSDA.pnts@data[sel,"TAXOUSDA.f"] = levs[j]
  }
}
TAXOUSDA.pnts$TAXOUSDA.f <- as.factor(TAXOUSDA.pnts$TAXOUSDA.f)
summary(TAXOUSDA.pnts$TAXOUSDA.f)
TAXOUSDA.pnts <- TAXOUSDA.pnts[!is.na(TAXOUSDA.pnts$TAXOUSDA.f),]
summary(as.factor(TAXOUSDA.pnts$SOURCEDB))
## Remove all spatial duplicates:
TAXOUSDA.pnts$LOC_ID <- as.factor(paste(TAXOUSDA.pnts@coords[,1], TAXOUSDA.pnts@coords[,2], sep="_"))
## 7495 duplicates!
TAXOUSDA.pnts <- TAXOUSDA.pnts[!duplicated(TAXOUSDA.pnts$LOC_ID),]
tax <- levels(TAXOUSDA.pnts$TAXOUSDA.f)
pal <- col.legend[match(tax, col.legend$Group),"COLOR"]
shape = "http://maps.google.com/mapfiles/kml/paddle/wht-blank.png"
kml(TAXOUSDA.pnts, subfolder.name="Observed TAXOUSDA", colour=TAXOUSDA.f, colour_scale=pal, shape=shape, size=.8, kmz=TRUE, labels=TAXOUSDA.pnts$TAXOUSDA.f)
str(TAXOUSDA.pnts@data)
save("TAXOUSDA.pnts", file="TAXOUSDA.pnts.rda")
## 30,270 points!

#load("equi7t3.rda")
source("extract.equi7t3.R")
ov <- extract.equi7t3(x=TAXOUSDA.pnts, y=c("C01GLC5","DEMMRG5","ES1MOD5","M01MOD4"), equi7t3=equi7t3, path="G:/SoilGrids250m/covs", cpus=4)
