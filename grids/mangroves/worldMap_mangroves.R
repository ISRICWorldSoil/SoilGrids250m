## Rasterize world map of Mangrove distribution
## http://data.unep-wcmc.org/datasets/4
## Tom.Hengl@isric.org

library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
gdal.dir <- shortPathName("C:/Program files/GDAL")
gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe")
gdal_rasterize <- paste0(gdal.dir, "/gdal_rasterize.exe")
saga_cmd <- shortPathName("C:\\Program Files\\SAGA-GIS\\saga_cmd.exe")
system(paste(saga_cmd, '--version'))
system(paste(gdal_rasterize, './DownloadPack-WCMC-010-MangroveUSGS2011-Ver1-3/WCMC-010-MangroveUSGS2011-ver1-3.shp', '-l WCMC-010-MangroveUSGS2011-ver1-3', '-te -180 -90 180 90', '-tr 0.002083333 0.002083333 -ot Byte', '-burn 100 MANGPR_250m.tif -a_nodata 0 -co \"COMPRESS=DEFLATE\"'))

#system(paste(saga_cmd, '-c=10 grid_gridding 0 -INPUT=\"./DownloadPack-WCMC-010-MangroveUSGS2011-Ver1-3/WCMC-010-MangroveUSGS2011-ver1-3.shp\" -FIELD=0 -OUTPUT=0 -GRID=MANGPR_250m.sgrd -GRID_TYPE=0 -TARGET_DEFINITION=0 -TARGET_USER_XMIN=\"-180\" -TARGET_USER_XMAX=\"180\" -TARGET_USER_YMIN=\"-90\" -TARGET_USER_YMAX=\"90\" -TARGET_USER_FITS=1 -TARGET_USER_SIZE=\"0.002083333\"'))
#system(paste(gdalwarp, 'MANGPR_250m.sdat  -tr 0.002083333 0.002083333 -ot \"Byte\" -dstnodata 255'))