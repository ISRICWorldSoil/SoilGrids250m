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

###################################
## Pseudo-points
###################################

## Desert / barerock soils (see "Encyclopedia of Soil Science" Vol.1 P520)
#deserts <- raster("../../masks/desertPR_sin.tif")
#barerock <- raster("../../masks/barerockPR_sin.tif")
#deserts.pnt <- sampleRandom(deserts, size=2000, sp=TRUE)
#plot(deserts.pnt) ## ca 400 points
#save(deserts.pnt, file="deserts.pnt.rda")
#barerock.pnt <- sampleRandom(barerock, size=30000, sp=TRUE)
#plot(barerock.pnt) ## ca 200 points 
#save(barerock.pnt, file="barerock.pnt.rda")
load("deserts.pnt.rda")
load("barerock.pnt.rda")

## Nicaragua points (http://hdl.handle.net/1902.1/20164)
#Nicaraqua.pol <- readOGR("../../../soilstorage/PolygonMaps/Nicaragua/taxonnic_wgs84.shp", "taxonnic_wgs84")
#Nicaraqua.sim <- spsample(Nicaraqua.pol, type="random", n=150)
#Nicaraqua.ov <- cbind(over(Nicaraqua.sim, Nicaraqua.pol["SUB_GRUPO"]), Nicaraqua.sim@coords)
#Nicaraqua.ov <- Nicaraqua.ov[!is.na(Nicaraqua.ov$SUB_GRUPO),]
#names(Nicaraqua.ov)[1] <- "TAXOUSDA"
#save(Nicaraqua.ov, file="Nicaraqua.ov.rda")
load("Nicaraqua.ov.rda")

## Circum-polar_regions points (http://bolin.su.se/data/ncscd/)
#polar.pol <- readOGR("../../../soilstorage/PolygonMaps/Circum-polar_regions/LAEA/NCSCD_Circumarctic_LAEA.shp", "NCSCD_Circumarctic_LAEA")
#turbels.pnt <- spsample(polar.pol[polar.pol$TURBEL_PCT>60,], type="random", n=300)
#orthels.pnt <- spsample(polar.pol[polar.pol$ORTHEL_PCT>60,], type="random", n=200)
#histels.pnt <- spsample(polar.pol[polar.pol$HISTEL_PCT>60,], type="random", n=100)
#save(turbels.pnt, file="turbels.pnt.rda")
#save(orthels.pnt, file="orthels.pnt.rda")
#save(histels.pnt, file="histels.pnt.rda")
load("turbels.pnt.rda")
load("orthels.pnt.rda")
load("histels.pnt.rda")

## Glacier areas from (http://www.glims.org/download/):
data(landmask20km)
load("glaciers.pnt.rda")
## suborders from the "NRCS Global Soil Regions Map" 
glaciers.pnt$suborder <- paste(over(y=landmask20km["suborder"], x=spTransform(glaciers.pnt, landmask20km@proj4string))$suborder)
sel.gelands = glaciers.pnt$suborder=="Ocean" & glaciers.pnt$SLPSRM3a <0.3 & glaciers.pnt$geog_area == "Iceland" | glaciers.pnt$geog_area == "Various (GlobGlacier)"
gelands.pnt <- as(glaciers.pnt[sel.gelands,], "SpatialPoints")
save(gelands.pnt, file="gelands.pnt.rda")
sel.gelepts = glaciers.pnt$suborder=="Gelepts" 
gelepts.pnt <- as(glaciers.pnt[sel.gelepts,], "SpatialPoints")
save(gelepts.pnt, file="gelepts.pnt.rda")
sel.gelods = glaciers.pnt$suborder=="Gelods"
gelods.pnt <- as(glaciers.pnt[sel.gelods,], "SpatialPoints")
save(gelods.pnt, file="gelods.pnt.rda")
sel.turbels = glaciers.pnt$suborder=="Turbels"|glaciers.pnt$suborder=="Ice" 
turbels2.pnt <- as(glaciers.pnt[sel.turbels,], "SpatialPoints")
save(turbels2.pnt, file="turbels2.pnt.rda")
sel.orthels = glaciers.pnt$suborder=="Orthels"|glaciers.pnt$suborder=="Ocean" & glaciers.pnt$SLPSRM3a >0.3
orthels2.pnt <- as(glaciers.pnt[sel.orthels,], "SpatialPoints")
save(orthels2.pnt, file="orthels2.pnt.rda")

## Indonesia points:
#ind.pol <- readOGR("../../../soilstorage/PolygonMaps/Indonesia/Kalimantan/tot_geo_shape.shp", "tot_geo_shape")
#summary(ind.pol$LS_CODE)
## attribute table:
#ind.code <- rbind.fill(lapply(paste0("../../../soilstorage/PolygonMaps/Indonesia/Kalimantan/facet6", 1:4, ".DBF"), function(x){read.dbf(x)$dbf}))
#summary(ind.code$L_SYS)
#ind.code <- plyr::rename(ind.code, replace=c("L_SYS"="LS_CODE"))
#ind.pol$TAXOUSDA <- join(ind.pol@data, ind.code, by="LS_CODE", type="left", match="first")$SOIL_NAME
#Indonesia.sim <- spsample(ind.pol, type="random", n=150)
#Indonesia.ov <- cbind(over(Indonesia.sim, ind.pol["TAXOUSDA"]), Indonesia.sim@coords)
#save(Indonesia.ov, file="Indonesia.ov.rda")
#x = summary(Indonesia.ov$TAXOUSDA)
#write.csv(x, file="cleanup_Indonesia.csv")
load("Indonesia.ov.rda")

n.lst <- c("Psamments", "Orthents", "Turbels", "Orthels", "Histels", "Gelands", "Gelepts", "Gelods", "Turbels", "Orthels") ## 
TAXOUSDA.sim <- lapply(list(as(deserts.pnt, "SpatialPoints"), as(barerock.pnt, "SpatialPoints"), turbels.pnt, orthels.pnt, histels.pnt, gelands.pnt, gelods.pnt, turbels2.pnt, orthels2.pnt), spTransform, CRS("+proj=longlat +datum=WGS84"))
TAXOUSDA.sim.df <- list(NULL)
for(j in 1:length(TAXOUSDA.sim)){
   TAXOUSDA.sim.df[[j]] <- cbind(as.data.frame(TAXOUSDA.sim[[j]]), TAXOUSDA=rep(n.lst[j], length(TAXOUSDA.sim[[j]])))
}
TAXOUSDA.sim.df[[length(n.lst)+1]] <- plyr::rbind.fill(Nicaraqua.ov, Indonesia.ov)
TAXOUSDA.sim <- plyr::rbind.fill(TAXOUSDA.sim.df)
str(TAXOUSDA.sim) ## 1406 pseudo points
TAXOUSDA.sim <- plyr::rename(TAXOUSDA.sim, c("x"="LONWGS84", "y"="LATWGS84"))
TAXOUSDA.sim$SOURCEID <- paste("SIM", 1:nrow(TAXOUSDA.sim), sep="_")
TAXOUSDA.sim$SOURCEDB = "Simulated"
#plot(TAXOUSDA.sim[,1:2])
## add simulated points to the list:
in.lst[[length(tax.lst)+1]] <- TAXOUSDA.sim

###################################
## Bind together and clean up names
###################################

all.pnts <- dplyr::rbind_all(in.lst)
str(all.pnts)
## 58,942
all.pnts <- as.data.frame(all.pnts)
coordinates(all.pnts) <- ~ LONWGS84+LATWGS84
proj4string(all.pnts) <- "+proj=longlat +datum=WGS84"
## Remove spatial duplicates (except for the NCSS data):
all.pnts$LOC_ID <- as.factor(paste(all.pnts@coords[,1], all.pnts@coords[,2], sep="_"))
summary(!duplicated(all.pnts$LOC_ID))
## >10,000 duplicate points
selP <- all.pnts$SOURCEDB=="NCSS"
summary(selP)
#selP.pnts <- all.pnts[selP,]
sel.pnts <- all.pnts[!selP,]
## remove duplicates, but keep all NCSS points
#selP.pnts <- selP.pnts[!duplicated(selP.pnts$LOC_ID),]
sel.pnts <- sel.pnts[!duplicated(sel.pnts$LOC_ID),]
all.pnts <- rbind(all.pnts[selP,], sel.pnts)
length(all.pnts)
## 54926
levels(as.factor(all.pnts$SOURCEDB))

## Targeted classes:
USDA_levs <- read.csv("USDA_levs.csv", na.strings = c("NA",""))
str(USDA_levs)
levs <- levels(USDA_levs$Group)[-which(levels(USDA_levs$Group) %in% c("Shifting Sand", "Rock", "Ocean", "Ice"))]
## 70 suborders
levsGG <- read.csv("TAXOUSDA_GreatGroups.csv")
#summary(levsGG$Great_Group)
strip_s <- function(x){ ifelse(substr(x, nchar(x), nchar(x))=="s", substr(x, 1, nchar(x)-1), x) }

## Suborder / Great_Group names --> we try to locate them in the raw names one by one
taxn.lst <- list(NULL)
## takes 2 mins:
for(j in 1:length(levs)){
  ## remove "s" if at the end of the class name:
  pat <- strip_s(levs[j]) 
  sel1 <- grep(pat, paste(all.pnts$TAXOUSDA), ignore.case=TRUE)
  ## also check in the full name - Subgroup name:
  ## Example: "Euic, pergelic Fluvaquentic Sapristel" is a "Histel", NOT a "Aquent"
  ## "Fine loamy,mixed,mesic,Xerofluventic Haplocambids" is NOT a "Fluvent"
  pat2 <- strip_s(paste(levsGG[which(levsGG$Suborder %in% levs[j]),"Great_Group"]))
  sel2 <- unlist(sapply(pat2, function(x){grep(x, paste(all.pnts$TAXNUSDA), ignore.case=TRUE)}))
  sel3 <- unlist(sapply(pat2, function(x){grep(x, paste(all.pnts$TAXOUSDA), ignore.case=TRUE)}))
  sel4 <- which(all.pnts$TAXOUSDA == levs[j])
  sel <- unique(c(sel1, sel2, sel3, sel4))
  ## bind together
  ## there can be multiple soil orders at the same location
  if(length(sel)>0){
    taxn.lst[[j]] <- data.frame(all.pnts[sel,])
    taxn.lst[[j]]$TAXOUSDA.f <- levs[j]
  }
}
TAXOUSDA.pnts <- do.call(rbind, taxn.lst)
TAXOUSDA.pnts$TAXOUSDA.f <- as.factor(TAXOUSDA.pnts$TAXOUSDA.f)
## examples:
TAXOUSDA.pnts[TAXOUSDA.pnts$SOURCEID=="S1995AK185018",]
TAXOUSDA.pnts[TAXOUSDA.pnts$SOURCEID=="IranSoil_54.48415_36.85611111",]
summary(TAXOUSDA.pnts$TAXOUSDA.f)
c.TAXOUSDA.f <- summary(TAXOUSDA.pnts$TAXOUSDA.f)
## Cryids and Gelepts <10 observations
write.csv(c.TAXOUSDA.f, "summary_TAXOUSDA.csv")

###################################
## Export points
###################################

length(TAXOUSDA.pnts$TAXOUSDA.f)
## FINAL NUMBER OF POINTS: 58,124 points
summary(as.factor(TAXOUSDA.pnts$SOURCEDB))
#    Can_FECD                                    CanSIS 
#         1086                                      6545 
#         ISIS                                      NCSS 
#          445                                     31468 
#    Russia_EGRPR                                USGS_Buell 
#          486                                      1060 
#         WISE                                    AfSPDB 
#         2233                                      1700 
#        CIFOR                                    eSOTER 
#          212                                       109 
#    RadamBrasil                                    Alaska 
#         6305                                       371 
#    Simulated                                  IranSoil 
#         2433                                      1314 
#    AU_NatSoil                                     Artic 
#         2076                                       275 
coordinates(TAXOUSDA.pnts) <- ~ LONWGS84+LATWGS84
proj4string(TAXOUSDA.pnts) <- "+proj=longlat +datum=WGS84"
unlink("TAXOUSDA.pnts.shp")
writeOGR(TAXOUSDA.pnts, "TAXOUSDA.pnts.shp", "TAXOUSDA.pnts", "ESRI Shapefile")
#plotKML(TAXOUSDA.pnts["TAXOUSDA.f"], file.name="TAXOUSDA_Feb_12_2016.kml", kmz=TRUE)
#str(TAXOUSDA.pnts@data)
#TAXOUSDA.pnts <- TAXOUSDA.pnts[c("SOURCEDB","SOURCEID","TAXOUSDA.f")]
save(TAXOUSDA.pnts, file="TAXOUSDA.pnts.rda")
