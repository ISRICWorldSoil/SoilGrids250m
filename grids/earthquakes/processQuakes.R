## Derive global density of earthquakes using:
## National Geophysical Data Center / World Data Service (NGDC/WDS): Significant Earthquake Database. National Geophysical Data Center, NOAA. doi:10.7289/V5TD9V7K
## USGS Earthquake Archives http://earthquake.usgs.gov/earthquakes/

library(sp)
library(rgdal)
library(raster)
library(spatstat)
library(maptools)
library(plyr)
gdalwarp = "/usr/local/bin/gdalwarp"
gdal_translate = "/usr/local/bin/gdal_translate"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")

## read quakes and focus on intensity:
quakes <- rbind.fill(lapply(list.files(pattern=".csv"), read.csv, stringsAsFactors=FALSE))
#plot(quakes[,c("longitude","latitude")])
#x <- read.delim("signif_quakes.txt")
x <- as.data.frame(readOGR("quakes.shp", "quakes"))
x <- plyr::rename(x, replace=c("coords.x1"="longitude", "coords.x2"="latitude", "YEAR"="time", "MAGF"="mag"))
sel = c("longitude","latitude","mag","time")
quakes <- rbind(quakes[,sel], x[,sel])
quakes$mag <- ifelse(is.na(quakes$mag), 6, quakes$mag)

GDALinfo("/data/MOD11A2/landMask1km_B.sdat")
unlink("landMask2km.tif")
system(paste0(gdalwarp, ' /data/MOD11A2/landMask1km_B.sdat landMask2km.tif -r \"near\" -tr 0.025 0.025 -co \"COMPRESS=DEFLATE\"'))
#system(paste0(gdalwarp, ' /data/MOD11A2/landMask1km_B.sdat landMask5km.tif -r \"near\" -tr 0.05 0.05 -co \"COMPRESS=DEFLATE\"'))
mask <- readGDAL("landMask2km.tif")
mask$total = 1
wowin <- as(mask["total"], "owin")
#wowin <- owin(c(-180,180), c(-90,90))
quakes.ppp <- ppp(quakes$longitude, quakes$latitude, marks=quakes$mag, window=wowin)
#plot(quakes.ppp)
## TAKES >30mins:
densMAG <- density.ppp(quakes.ppp, sigma=0.5, weights=quakes.ppp$marks)
dens.MAG = as(densMAG, "SpatialGridDataFrame")
writeGDAL(dens.MAG, "dens.MAG.tif", type="Int16", options="COMPRESS=DEFLATE")
