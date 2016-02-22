## Crowdsourced data UKSO
## Acidity and soil texture

#library(gdalUtils)
library(rgdal)
gdal.dir = shortPathName("C:\\Program Files\\GDAL")
ogr2ogr = paste(gdal.dir, "ogr2ogr.exe", sep="\\")
gdalinfo = paste(gdal.dir, "gdalinfo.exe", sep="\\")
gdal_translate = paste(gdal.dir, "gdal_translate", sep="\\")

system(paste0(gdalinfo, ' WMS:https://map.bgs.ac.uk/arcgis/services/UKSO/UKSO_Crowdsourced/MapServer/WmsServer?'))
## download data from WMS:
#system(paste0(ogr2ogr, ' -f \"ESRI Shapefile\UKSO_texture.shp WMS:\"https://map.bgs.ac.uk/arcgis/services/UKSO/UKSO_Crowdsourced/MapServer/WmsServer?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=Crowdsourced.Soil.Texture&SRS=EPSG:4326&BBOX=-180.000000,-90.000000,180.000000,90.000000"'))
system(paste0(gdal_translate, ' "https://map.bgs.ac.uk/arcgis/services/UKSO/UKSO_Crowdsourced/MapServer/WmsServer?f=json" wms_ukso.xml -of WMS'))
