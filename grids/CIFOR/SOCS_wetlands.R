## Assessment and modeling of Soil Organic Carbon Stocks for Wetlands and Peatlands in the tropics
## Tom.Hengl@isric.org

library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
library(plotKML)
library(psych)
library(scales)
library(R.utils)
library(plotKML)
library(GSIF)
plotKML.env(convert="convert", show.env=FALSE)

if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate = "/usr/bin/gdal_translate"
  gdalwarp = "/usr/bin/gdalwarp"
}

## List of property maps:
fao.lst <- c("Sapric.Histosols", "Hemic.Histosols", "Fibric.Histosols", "Cryic.Histosols", "Histic.Albeluvisols")
usda.lst <- c("Saprists", "Hemists", "Folists", "Fibrists")

## resample to 1 km resolution Tropics only:
system(paste0(gdalwarp, ' /data/GEOG/HISTPR_1km_ll.tif TROP_HISTPR_1km_ll.tif -te -180 -36 180 36 -co \"COMPRESS=DEFLATE\"'))
for(i in fao.lst){ system(paste0(gdalwarp, ' /data/GEOG/TAXNWRB_', i, '_1km_ll.tif TROP_', i, '_1km_ll.tif -te -180 -36 180 36 -co \"COMPRESS=DEFLATE\"')) }

## Resample in parallel:
orc.lst <- paste0("OCSTHA_M_sd", 1:6)
bld.lst <- paste0("BLDFIE_M_sl", 1:7)
soc.lst <- paste0("ORCDRC_M_sl", 1:7)
t.lst <- c(orc.lst, bld.lst, soc.lst)
sfInit(parallel=TRUE, cpus=length(t.lst))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("gdalwarp","t.lst")
x <- sfClusterApplyLB(t.lst, function(x){ try( if(!file.exists(paste0('TROP_', x, '_1km_ll.tif'))){ system(paste0(gdalwarp, ' /data/GEOG/', x, '_1km_ll.tif TROP_', x, '_1km_ll.tif -te -180 -36 180 36 -co \"COMPRESS=DEFLATE\"')) } ) })
sfStop()
file.copy(from="TROP_OCSTHA_M_sd6_1km_ll.tif", to="TROP_OCSTHA_2m_1km_ll.tif")
## sum up OCS values for 0-1 m
sD <- raster::stack(paste0('TROP_OCSTHA_M_sd', 1:5, '_1km_ll.tif'))
sumf <- function(x){calc(x, sum, na.rm=TRUE)}
## run in parallel:
beginCluster()
r1 <- clusterR(sD, fun=sumf, filename="TROP_OCSTHA_1m_1km_ll.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
endCluster()

#for(i in orc.lst){ system(paste0(gdalwarp, ' /data/GEOG/', i, '_1km_ll.tif TROP_', i, '_1km_ll.tif -te -180 -36 180 36 -co \"COMPRESS=DEFLATE\"')) }

## plot HWSD and SoilGrids estimated SOCS next to each other:
rob = "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
system(paste0(gdalwarp, ' TROP_OCSTHA_1m_1km_ll.tif OCSTHA_1m_10km.tif -r \"average\" -tr 0.1 0.1 -co \"COMPRESS=DEFLATE\" -dstnodata -9999'))
system(paste0(gdalwarp, ' TROP_OCSTHA_2m_1km_ll.tif OCSTHA_2m_10km.tif -r \"average\" -tr 0.1 0.1 -co \"COMPRESS=DEFLATE\" -dstnodata -9999'))
system(paste0(gdalwarp, ' TROP_HISTPR_1km_ll.tif HISTPR_10km.tif -r \"average\" -tr 0.1 0.1 -co \"COMPRESS=DEFLATE\" -dstnodata 255'))

## Soil points:
load("/data/models/SPROPS/ovA.rda")
TROP_xy <- ovA[ovA$LATWGS84>-36&ovA$LATWGS84<36,c("SOURCEID","LONWGS84","LATWGS84","SOURCEDB","UHDICM","LHDICM","HZDTXT","ORCDRC","BLD","CRFVOL")]
TROP_xy <- TROP_xy[!is.na(TROP_xy$LATWGS84)&!is.na(TROP_xy$ORCDRC),]
unlink("TROP_soil_profiles.csv.gz")
write.csv(TROP_xy, file="TROP_soil_profiles.csv")
gzip("TROP_soil_profiles.csv")
summary(as.factor(TROP_xy$SOURCEDB))
plot(TROP_xy$LONWGS84, TROP_xy$LATWGS84, pch="+", col="red")
rm(ovA)

load("/data/models/TAXOUSDA/ov.TAXOUSDA.rda")
TROP_usda <- ov[ov$LATWGS84>-36&ov$LATWGS84<36,c("SOURCEID","LOC_ID","LATWGS84","SOURCEDB","TAXOUSDA.f")]
TROP_usda$LONWGS84 <- as.numeric(sapply(paste(TROP_usda$LOC_ID), function(i){strsplit(i, "_")[[1]][1]}))
TROP_usda$HISTPR <- 0
TROP_usda$HISTPR[grep(TROP_usda$TAXOUSDA.f, pattern="ist", ignore.case=TRUE)] <- 1
summary(as.factor(TROP_usda$HISTPR))
TROP_usda <- TROP_usda[!is.na(TROP_usda$LATWGS84)&!is.na(TROP_usda$TAXOUSDA.f),]
unlink("TROP_soil_types_USDA.csv.gz")
write.csv(TROP_usda[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXOUSDA.f","HISTPR")], file="TROP_soil_types_USDA.csv")
gzip("TROP_soil_types_USDA.csv")
summary(TROP_usda$SOURCEDB)
points(TROP_usda$LONWGS84, TROP_usda$LATWGS84, pch="+")

load("/data/models/TAXNWRB/ov.TAXNWRB.rda")
TROP_wrb <- ov[ov$LATWGS84>-36&ov$LATWGS84<36,c("SOURCEID","LOC_ID","LATWGS84","SOURCEDB","TAXNWRB.f")]
TROP_wrb$LONWGS84 <- as.numeric(sapply(paste(TROP_wrb$LOC_ID), function(i){strsplit(i, "_")[[1]][1]}))
TROP_wrb$HISTPR <- 0
TROP_wrb$HISTPR[grep(TROP_wrb$TAXNWRB.f, pattern="hist", ignore.case=TRUE)] <- 1
summary(as.factor(TROP_wrb$HISTPR))
TROP_wrb <- TROP_wrb[!is.na(TROP_wrb$LATWGS84)&!is.na(TROP_wrb$TAXNWRB.f),]
unlink("TROP_soil_types_WRB.csv.gz")
write.csv(TROP_wrb[,c("SOURCEID","SOURCEDB","LONWGS84","LATWGS84","TAXNWRB.f","HISTPR")], file="TROP_soil_types_WRB.csv")
gzip("TROP_soil_types_WRB.csv")
rm(ov)

## Wetlands / land cover classes:
glc.lst <- paste0("/data/GlobCover30/", c(paste0("L0",1:9,"GLC3a.tif"), "L10GLC3a.tif", "LMKGLC3a.tif"))
sfInit(parallel=TRUE, cpus=length(glc.lst))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("gdalwarp","glc.lst")
x <- sfClusterApplyLB(glc.lst, function(x){ try( system(paste0(gdalwarp, ' ', x, ' ', gsub("3a", "_10km", basename(x)),' -r \"average\" -tr 0.1 0.1 -te -180 -36 180 36 -co \"COMPRESS=DEFLATE\"')) ) })
sfStop()

## Upper and lower limits for OCS at 10 km:
tif.lst <- paste0("/data/GEOG/", c(paste0("ORCDRC_M_sl", 1:7, "_1km_ll.tif"), paste0("BLDFIE_M_sl", 1:7, "_1km_ll.tif"), paste0("CRFVOL_M_sl", 1:7, "_1km_ll.tif")))
sfInit(parallel=TRUE, cpus=length(tif.lst))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("gdalwarp","tif.lst")
x <- sfClusterApplyLB(tif.lst, function(x){ try( system(paste0(gdalwarp, ' ', x, ' ', gsub("1km", "10km", basename(x)),' -r \"average\" -tr 0.1 0.1 -te -180 -36 180 36 -co \"COMPRESS=DEFLATE\"')) ) })
sfStop()

## horizon thickness:
ds <- get("stsize", envir = GSIF.opts)
ds <- rowMeans(data.frame(c(NA,ds),c(ds,NA)), na.rm=TRUE)
g10km <- raster::stack(list.files(pattern="10km"))
names(g10km)
g10km <- as(g10km, "SpatialGridDataFrame")
## Weighted average based on the thickness of horizon:
g10km$BLDFIE_M_1m_10km_ll <- rowSums(g10km@data[,paste0("BLDFIE_M_sl", 1:6, "_10km_ll")] * data.frame(lapply(ds[-7], rep, length=nrow(g10km))), na.rm=TRUE) / sum(ds[-7])
g10km$ORCDRC_M_1m_10km_ll <- rowSums(g10km@data[,paste0("ORCDRC_M_sl", 1:6, "_10km_ll")] * data.frame(lapply(ds[-7], rep, length=nrow(g10km))), na.rm=TRUE) / sum(ds[-7])
g10km$CRFVOL_M_1m_10km_ll <- rowSums(g10km@data[,paste0("CRFVOL_M_sl", 1:6, "_10km_ll")] * data.frame(lapply(ds[-7], rep, length=nrow(g10km))), na.rm=TRUE) / sum(ds[-7])
g10km$ORCDRC_M_2m_10km_ll <- rowSums(g10km@data[,c("ORCDRC_M_sl6_10km_ll","ORCDRC_M_sl7_10km_ll")], na.rm=TRUE)
g10km$BLDFIE_M_2m_10km_ll <- rowSums(g10km@data[,c("BLDFIE_M_sl6_10km_ll","BLDFIE_M_sl7_10km_ll")], na.rm=TRUE)
g10km$CRFVOL_M_2m_10km_ll <- rowSums(g10km@data[,c("CRFVOL_M_sl6_10km_ll","CRFVOL_M_sl7_10km_ll")], na.rm=TRUE)
plot(raster(g10km["BLD_M_sd1_10km_ll"]), col=SAGA_pal[[1]])
g10km <- as(g10km, "SpatialPixelsDataFrame")
g10km$HISTPR_10km <- readGDAL("HISTPR_10km.tif")$band1[g10km@grid.index]
summary(!is.na(g10km$LMKGLC_10km))
summary(!is.na(g10km$OCSTHA_2m_10km))
g10km <- g10km[!is.na(g10km$LMKGLC_10km),]
plot(raster(g10km["BLDFIE_M_sl1_10km_ll"]), col=SAGA_pal[[1]])
plot(raster(g10km["HISTPR_10km"]), col=SAGA_pal[[1]])
plot(raster(g10km["BLDFIE_M_1m_10km_ll"]), col=SAGA_pal[[1]])
plot(log1p(raster(g10km["OCSTHA_1m_10km"])), col=SAGA_pal[[1]])
plot(log1p(raster(g10km["OCSTHA_2m_10km"])), col=SAGA_pal[[1]])
plot(log1p(raster(g10km["ORCDRC_M_1m_10km_ll"])), col=SAGA_pal[[1]])

## UPPER / LOWER UNCERTAINTY ESTIMATES:
OCS_1m <- GSIF::OCSKGM(ORCDRC=g10km$ORCDRC_M_1m_10km_ll, BLD=g10km$BLDFIE_M_1m_10km_ll, CRFVOL=g10km$CRFVOL_M_1m_10km_ll, HSIZE=100, ORCDRC.sd=20, BLD.sd=170) ## (expm1(log1p(g10km$ORCDRC_M_1m_10km_ll)+0.6)-expm1(log1p(g10km$ORCDRC_M_1m_10km_ll)-0.6))/2
OCS_2m <- GSIF::OCSKGM(ORCDRC=g10km$ORCDRC_M_2m_10km_ll, BLD=g10km$BLDFIE_M_2m_10km_ll, CRFVOL=g10km$CRFVOL_M_2m_10km_ll, HSIZE=100, ORCDRC.sd=20, BLD.sd=170)
g10km$OCSTHA_1m_10km_UPPER <- g10km$OCSTHA_1m_10km + attr(OCS_1m, "measurementError")*10
g10km$OCSTHA_1m_10km_LOWER <- g10km$OCSTHA_1m_10km - attr(OCS_1m, "measurementError")*10
g10km$OCSTHA_1m_10km_LOWER <- ifelse(g10km$OCSTHA_1m_10km_LOWER<0, 0, g10km$OCSTHA_1m_10km_LOWER)
g10km$OCSTHA_2m_10km_UPPER <- g10km$OCSTHA_2m_10km + attr(OCS_2m, "measurementError")*10
g10km$OCSTHA_2m_10km_LOWER <- g10km$OCSTHA_2m_10km - attr(OCS_2m, "measurementError")*10
g10km$OCSTHA_2m_10km_LOWER <- ifelse(g10km$OCSTHA_2m_10km_LOWER<0, 0, g10km$OCSTHA_2m_10km_LOWER)
## test it:
g10km@data[119000,c("OCSTHA_1m_10km_LOWER","OCSTHA_1m_10km","OCSTHA_1m_10km_UPPER")]
g10km@data[which(g10km$OCSTHA_1m_10km>2200)[1],c("OCSTHA_1m_10km_LOWER","OCSTHA_1m_10km","OCSTHA_1m_10km_UPPER")]

g10km.pol <- grid2poly(g10km["LMKGLC_10km"]) ## Takes ca 10 mins!
library(geosphere)
## Calculate area in ha for each pixel in latlon (run in parallel):
getArea <- function(x,pol){geosphere::areaPolygon(as(pol[x,], "SpatialPolygons"))/1e4}
getArea(2,g10km.pol)
AREA <- unlist(parallel::mclapply( 1:length(g10km.pol), getArea, pol=g10km.pol, mc.cores=48))
#g10km$AREA <- sapply(1:length(g10km.pol), function(x){geosphere::areaPolygon(as(g10km.pol[x,], "SpatialPolygons"))})/1e4
g10km$AREA <- AREA
summary(g10km$AREA)

g10km.df = as.data.frame(g10km[,c("OCSTHA_1m_10km","OCSTHA_2m_10km","HISTPR_10km","L05GLC_10km","BLDFIE_M_1m_10km_ll","BLDFIE_M_2m_10km_ll","ORCDRC_M_1m_10km_ll","ORCDRC_M_2m_10km_ll","OCSTHA_1m_10km_LOWER","OCSTHA_1m_10km_UPPER","OCSTHA_2m_10km_LOWER","OCSTHA_2m_10km_UPPER","AREA")])
unlink("TROP_grid10km.csv.gz")
write.csv(g10km.df, file="TROP_grid10km.csv")
gzip("TROP_grid10km.csv")
## Kalimantan:
kal = which(g10km.df$s1>110.6 & g10km.df$s1<110.8 & g10km.df$s2> -2.9 & g10km.df$s2 < -2.7)
g10km.df[kal,c("OCSTHA_1m_10km_LOWER","OCSTHA_1m_10km","OCSTHA_1m_10km_UPPER")]

## plot in Google Earth:
setwd("/data/CIFOR")
#source("plotKML.GDALobj.R")
#source("legend.bar.R")
#r1 = raster("TROP_OCSTHA_1m_1km_ll.tif")
#r2 = raster("TROP_OCSTHA_2m_1km_ll.tif")
#r3 = raster("HISTPR_1km_ll.tif")
#hist(log1p(sampleRandom(r1, 1e3)))
#setwd("/data/CIFOR/OCSTHA_2m")
#beginCluster()
#r <- clusterR(r2, fun=log1p, filename="TROP_OCSTHA_2m_1km_ll_log1p.tif", options=c("COMPRESS=DEFLATE"))
#endCluster()
#writeRaster(log1p(r), "TROP_OCSTHA_2m_1km_ll_log1p.tif")
#obj = GDALinfo("TROP_OCSTHA_2m_1km_ll_log1p.tif")
#plotKML.GDALobj(obj, file.name="TROP_log_OCSTHA_2m_1km.kml", block.x=5, z.lim=c(0,7.2), colour_scale = SAGA_pal[[1]], CRS="+proj=longlat +datum=WGS84", plot.legend=FALSE) # z.lim=c(0,1200)

#setwd("/data/CIFOR/OCSTHA_1m")
#beginCluster()
#r <- clusterR(r1, fun=log1p, filename="TROP_OCSTHA_1m_1km_ll_log1p.tif", options=c("COMPRESS=DEFLATE"))
#endCluster()
#obj = GDALinfo("TROP_OCSTHA_1m_1km_ll_log1p.tif")
#plotKML.GDALobj(obj, file.name="TROP_log_OCSTHA_1m_1km.kml", block.x=5, z.lim=c(0,7.2), colour_scale = SAGA_pal[[1]], CRS="+proj=longlat +datum=WGS84", plot.legend=FALSE)

#setwd("/data/CIFOR/HISTPROB")
#obj = GDALinfo("HISTPR_1km_ll.tif")
#plotKML.GDALobj(obj, file.name="TROP_HISTPR_1km.kml", block.x=5, z.lim=c(0,40), colour_scale = SAGA_pal[["SG_COLORS_YELLOW_BLUE"]], CRS="+proj=longlat +datum=WGS84", plot.legend=FALSE)
#setwd("/data/CIFOR")

## -tr 9000 9000 -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"', rob, '\" -co \"COMPRESS=DEFLATE\" -dstnodata -9999')) 
te <- as.vector(extent(raster("OCSTHA_1m_10km.tif")))
system(paste0(gdalwarp, ' X:/HWSD/HWSDa_OC_Dens_Top_5min.rst OCSTHA_1m_HWSDa.tif -tr 0.1 0.1 -co \"COMPRESS=DEFLATE\" -te ', paste(te[c(1,3,2,4)], collapse=" "))) 
##  -tr 9000 9000 -s_srs \"+proj=longlat +datum=WGS84\" -te ', paste(te[c(1,3,2,4)], collapse=" "),' -t_srs \"', rob, '\" -co \"COMPRESS=DEFLATE\"'))
system(paste0(gdalwarp, ' X:/HWSD/HWSDa_OC_Dens_Sub_5min.rst OCSTHA_1m_HWSDb.tif -tr 0.1 0.1 -co \"COMPRESS=DEFLATE\" -te ', paste(te[c(1,3,2,4)], collapse=" ")))

trop <- stack(c("OCSTHA_1m_10km.tif","OCSTHA_1m_HWSDa.tif","OCSTHA_1m_HWSDb.tif"))
trop <- as(as(trop, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
trop$SOCS_old <- trop$OCSTHA_1m_HWSDa+trop$OCSTHA_1m_HWSDb
names(trop)
rn = c(10,1200) ## quantile(c(trop$OCSTHA_1m_9km, trop$SOCS_old), c(.01, .99), na.rm=TRUE)
rx = rev(as.character(round(c(round(rn[1], 0), NA, round(mean(rn), 0), NA, round(rn[2], 0)), 2)))
trop$SOCS_oldf <- ifelse(trop$SOCS_old<rn[1], rn[1], ifelse(trop$SOCS_old>rn[2], rn[2], trop$SOCS_old))
trop$SOCS_SGf <- ifelse(trop$OCSTHA_1m_10km<rn[1], rn[1], ifelse(trop$OCSTHA_1m_10km>rn[2], rn[2], trop$OCSTHA_1m_10km))

require(maptools)
require(maps)
library(rgeos)
country <- map('world', plot=FALSE, fill=TRUE)
IDs <- sapply(strsplit(country$names, ":"), function(x) x[1])
country = as(map2SpatialPolygons(country, IDs=IDs), "SpatialLines")
b_poly <- as(extent(c(-180,180,-36,36)), "SpatialPolygons")
country = gIntersection(country, b_poly, byid = T)
proj4string(country) = "+proj=longlat +datum=WGS84"
#country <- spTransform(country, CRS(rob))

png(file="Fig_SOCS_comparison.png", res=100, width=1200, height=1200*2*(36+36)/180/2)
#spplot(trop, col.regions=SAGA_pal[[1]])
par(mfrow=c(2,1))
par(mai=c(0,0,0,0), oma=c(0,0,0,0),xaxs='i', yaxs='i')
image(log1p(raster(trop["SOCS_oldf"])), col=SAGA_pal[[1]], zlim=log1p(rn), main="", axes=FALSE, xlab="", ylab="") # , cex.lab=.7, cex.axis=.7
lines(country)
legend("left", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.8)
image(log1p(raster(trop["SOCS_SGf"])), col=SAGA_pal[[1]], zlim=log1p(rn), main="", axes=FALSE, xlab="", ylab="")
lines(country)
legend("left", rx, fill=rev(SAGA_pal[[1]][c(1,5,10,15,20)]), horiz=FALSE, cex=.8)
dev.off()

save.image()

## Scatter plot histograms:
df.s <- trop@data[sample.int(length(trop),20000),]
with(df.s, scatter.hist(SOCS_old,OCSTHA_1m_10km, xlab="HWSD", ylab="SoilGrids", pch=19, col=alpha("lightblue", 0.6), cex=1.5))
