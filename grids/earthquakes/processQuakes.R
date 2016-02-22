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
t_srs <- "+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

## read quakes and focus on intensity:
quakes <- rbind.fill(lapply(list.files(pattern=".csv"), read.csv, stringsAsFactors=FALSE))
#plot(quakes[,c("longitude","latitude")])
x <- as.data.frame(readOGR("quakes.shp", "quakes"))
x <- plyr::rename(x, replace=c("coords.x1"="longitude", "coords.x2"="latitude", "YEAR"="time", "MAGF"="mag"))
sel = c("longitude","latitude","mag","time")
quakes <- rbind(quakes[,sel], x[,sel])
quakes$mag <- ifelse(is.na(quakes$mag), 6, quakes$mag)
quakes$ID <- paste(quakes$longitude, quakes$latitude, quakes$mag, quakes$time, sep="_")
quakes <- quakes[!duplicated(quakes$ID),]
str(quakes) ## 84395 quakes 5+
save(quakes, file="quakes.rda")
coordinates(quakes) <- ~ longitude + latitude
proj4string(quakes) = "+proj=longlat +ellps=WGS84"
## volcanoes (http://www.ngdc.noaa.gov/hazard/volcano.shtml):
volcanoes <- read.csv("volcanoes.csv")
volcanoes <- volcanoes[!is.na(volcanoes$Latitude),]
str(volcanoes)
save(volcanoes, file="volcanoes.rda")
coordinates(volcanoes) <- ~ Longitude + Latitude
proj4string(volcanoes) = "+proj=longlat +ellps=WGS84"

GDALinfo("/data/MOD11A2/landMask1km_B.sdat")
unlink("landMask5km.tif")
system(paste0(gdalwarp, ' /data/MOD11A2/landMask1km_B.sdat landMask5km.tif -r \"near\" -tr 0.05 0.05 -te -180 -70 180 84 -co \"COMPRESS=DEFLATE\"'))
## reproject to TM projection to minimize distortions:
GDALinfo("landMask5km.tif")
unlink("landMask5km_tm.tif")
system(paste0(gdalwarp, ' landMask5km.tif landMask5km_tm.tif -r \"near\" -tr 5000 5000 -t_srs \"', t_srs, '\" -co \"COMPRESS=DEFLATE\"'))
mask <- readGDAL("landMask5km_tm.tif")
mask$total = 1
wowin <- as(mask["total"], "owin")
quakes.xy <- spTransform(quakes, CRS(proj4string(mask)))
unlink("quakes_tm.shp")
writeOGR(quakes.xy, "quakes_tm.shp", "quakes_tm", "ESRI Shapefile")
#summary(wowin[["m"]][1,])
#summary(wowin[["m"]][,1])
#wowin <- owin(c(-180,180), c(-90,90))
quakes.ppp <- ppp(quakes.xy@coords[,1], quakes.xy@coords[,2], marks=quakes.xy$mag, window=wowin)
## 71 points were rejected as lying outside the specified window
#plot(quakes.ppp)
## TAKES >1hr:
densMAG <- density.ppp(quakes.ppp, sigma=30000, weights=quakes.ppp$marks)
dens.MAG = as(densMAG, "SpatialGridDataFrame")
dens.MAG$vf = dens.MAG$v*1e9
#plot(raster(dens.MAG["vf"]), col=rev(bpy.colors(30)))
unlink("dens.MAG.tif")
writeGDAL(dens.MAG["vf"], "dens.MAG.tif", type="Int16", options="COMPRESS=DEFLATE")
## This one takes only few minutes, but shows some artifacts (circles)!
#system(paste0('/usr/local/bin/saga_cmd -c=48 grid_gridding 6 -POINTS=\"quakes_tm.shp\" -POPULATION=\"mag\" -RADIUS=60000 -KERNEL=1 -TARGET_DEFINITION=0 -TARGET_USER_XMIN=', mask@bbox["x","min"],' -TARGET_USER_XMAX=', mask@bbox["x","max"],' -TARGET_USER_YMIN=', mask@bbox["y","min"],' -TARGET_USER_YMAX=', mask@bbox["y","max"],' -TARGET_USER_SIZE=5000 -TARGET_USER_FITS=0 -TARGET_OUT_GRID=\"dens.MAG.sgrd\"'))
#system("7za a dens.MAG.7z dens.MAG.*")
#plot(raster("dens.MAG.sdat"))
#system(paste0(gdalwarp, ' dens.MAG.sdat dens.MAG_1km.tif -ot \"Int16\" -co \"COMPRESS=DEFLATE\" -dstnodata \"-32767\" -tr 1000 1000 -r \"cubicspline\"'))

volcanoes.xy <- spTransform(volcanoes, CRS(proj4string(mask)))
volcanoes.xy$marks = 1
unlink("volcanoes_tm.shp")
writeOGR(volcanoes.xy, "volcanoes_tm.shp", "volcanoes_tm", "ESRI Shapefile")
volcanoes.ppp <- ppp(volcanoes.xy@coords[,1], volcanoes.xy@coords[,2], marks=volcanoes.xy$mag, window=wowin)
## 16 points were rejected as lying outside the specified window
#plot(quakes.ppp)
## TAKES >1hr:
dens.volc <- density.ppp(volcanoes.ppp, sigma=60000)
dens.volc = as(dens.volc, "SpatialGridDataFrame")
dens.volc$vf = dens.volc$v*1e12
plot(raster(dens.volc["vf"]), col=rev(bpy.colors(30)))
unlink("volcanoes5km.tif")
writeGDAL(dens.volc["vf"], "volcanoes5km.tif", type="Int16", options="COMPRESS=DEFLATE")
#system(paste0('/usr/local/bin/saga_cmd -c=48 grid_gridding 6 -POINTS=\"volcanoes_tm.shp\" -POPULATION=\"marks\" -RADIUS=60000 -KERNEL=1 -TARGET_DEFINITION=0 -TARGET_USER_XMIN=', mask@bbox["x","min"],' -TARGET_USER_XMAX=', mask@bbox["x","max"],' -TARGET_USER_YMIN=', mask@bbox["y","min"],' -TARGET_USER_YMAX=', mask@bbox["y","max"],' -TARGET_USER_SIZE=5000 -TARGET_USER_FITS=0 -TARGET_OUT_GRID=\"volcanoes5km.sgrd\"'))
#system("7za a volcanoes5km.7z volcanoes5km.*")
#plot(raster("volcanoes5km.sdat"))
#points(volcanoes.xy)
