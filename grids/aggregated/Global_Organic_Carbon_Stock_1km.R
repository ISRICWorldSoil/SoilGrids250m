## Organic carbon stock 0-1 m and 1-2 m in MODIS projection system:
## Tom.Hengl@isric.org

library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
library(R.utils)
library(parallel)
library(doParallel)
library(foreign)
library(doSNOW)
library(doMC)

if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program Files/QGIS 2.16.0/bin")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe")
} else {
  gdal_translate = "/usr/bin/gdal_translate"
  gdalwarp = "/usr/bin/gdalwarp"
}

## MODIS land cover image at 1 km resolution
system(paste0(gdalwarp, ' /data/MCD12Q1/LandCover_2013001_L1_500m.tif /data/MCD12Q1/LandCover_2013001_L1_1km.tif -tr 1000 1000 -r \"near\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000'))
r <- raster("/data/MCD12Q1/LandCover_2013001_L1_1km.tif")
extent(r)
ncols = ncol(r)
nrows = nrow(r)
xllcorner = extent(r)[1]
yllcorner = extent(r)[3]
xurcorner = extent(r)[2]
yurcorner = extent(r)[4]
cellsize = res(r)[1]
NODATA_value = -32768

#system("7za e modis_sinusoidal.zip")
#world_sin = readOGR("countries_sinusoidal_world.shp", "countries_sinusoidal_world")
#proj4string(world_sin)
## GAUL 2015 data set
system("7za e g2015_2014_0.zip")

countries = read.dbf("g2015_2014_0.dbf")
str(countries$dbf$ADM0_NAME)
countries$dbf$ADM0_INT = as.integer(countries$dbf$ADM0_NAME)
write.dbf(countries, "g2015_2014_0.dbf")
countries.sum = data.frame(NAMES=levels(countries$dbf$ADM0_NAME), Value=1:length(levels(countries$dbf$ADM0_NAME)))
write.csv(countries.sum, "countries_legend.csv")
system(paste('ogr2ogr -t_srs \"', proj4string(r), '\" g2015_2014_0_a.shp g2015_2014_0.shp'))
system(paste0('/usr/local/bin/saga_cmd -c=48 grid_gridding 0 -INPUT \"g2015_2014_0_a.shp\" -FIELD \"ADM0_INT\" -GRID \"countries_100m.sgrd\" -GRID_TYPE 0 -TARGET_DEFINITION 0 -TARGET_USER_SIZE ', cellsize, ' -TARGET_USER_XMIN ', xllcorner+cellsize/2,' -TARGET_USER_XMAX ', xurcorner-cellsize/2, ' -TARGET_USER_YMIN ', yllcorner+cellsize/2,' -TARGET_USER_YMAX ', yurcorner-cellsize/2))
unlink("GAUL_COUNTRIES_1km")
system(paste(gdal_translate, 'countries_100m.sdat GAUL_COUNTRIES_1km.tif -ot \"Int16\" -co \"COMPRESS=DEFLATE\"'))


## Organic carbon stocks
orc.lst <- paste0("OCSTHA_M_sd", 1:6, "_1km_ll.tif")

sfInit(parallel=TRUE, cpus=length(orc.lst))
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("gdalwarp","orc.lst","r","cellsize")
x <- sfClusterApplyLB(orc.lst, function(x){ try( if(!file.exists(gsub("_ll", "_sin", x))){ system(paste0(gdalwarp, ' /data/GEOG/', x, ' ', gsub("_ll", "_sin", x), ' -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "))) } ) })
sfStop()

## sum up OCS values for 0-1 m
sD <- raster::stack(paste0('OCSTHA_M_sd', 1:5, '_1km_sin.tif'))
sumf <- function(x){calc(x, sum, na.rm=TRUE)}
## run in parallel:
beginCluster()
r1 <- clusterR(sD, fun=sumf, filename="OCSTHA_1m_1km_sin.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
endCluster()
file.copy(from="OCSTHA_M_sd6_1km_sin.tif", to="OCSTHA_2m_1km_sin.tif")
unlink("OCSTHA_M_sd6_1km_sin.tif")

## DERIVE Summary stats per country:
g1km = raster::stack(c("OCSTHA_1m_1km_sin.tif","OCSTHA_2m_1km_sin.tif","GAUL_COUNTRIES_1km.tif"))
g1km = as(g1km, "SpatialGridDataFrame")
sel.na = !is.na(g1km$GAUL_COUNTRIES_1km)
summary(sel.na)
g1km$Country_NAME = plyr::join(data.frame(Value=g1km$GAUL_COUNTRIES_1km), countries.sum, type="left")$NAMES
str(g1km@data)
summary(g1km$Country_NAME)
save.image()
gc(); gc()

registerDoMC(36)
SOC_agg <- ddply(g1km@data[sel.na,], .(Country_NAME), summarize, Total_OCS_1m_Pg=round(sum(OCSTHA_1m_1km_sin * 1e8,na.rm = TRUE)/1e15,1), Total_OCS_2m_Pg=round(sum(OCSTHA_2m_1km_sin * 1e8,na.rm = TRUE)/1e15,1), Total_OCS_1m_N=sum(!is.na(Country_NAME)), .parallel = TRUE)
closeAllConnections(); gc()
str(SOC_agg)
write.csv(SOC_agg, "Summary_OCS_per_country_1km.csv")

