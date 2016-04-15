## Rasterize world map of Mangrove distribution
## http://data.unep-wcmc.org/datasets/4
## Tom.Hengl@isric.org

library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
library(R.utils)
library(plotKML)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
gdal_rasterize <- "/usr/local/bin/gdal_rasterize"

unzip("WCMC-010-MangroveUSGS2011-Ver1-3.zip")
## CAN NOT BE PARALELLIZED:
#system(paste(gdal_rasterize, './DownloadPack-WCMC-010-MangroveUSGS2011-Ver1-3/WCMC-010-MangroveUSGS2011-ver1-3.shp', '-l WCMC-010-MangroveUSGS2011-ver1-3', '-te -180 -90 180 90', '-tr 0.002083333 0.002083333 -ot Byte', '-burn 100 MNGUSG_250m.tif -a_nodata 0 -co \"COMPRESS=DEFLATE\"'))

#system(paste('/usr/local/bin/saga_cmd -c=40 grid_gridding 0 -INPUT=\"./DownloadPack-WCMC-010-MangroveUSGS2011-Ver1-3/WCMC-010-MangroveUSGS2011-ver1-3.shp\" -FIELD=0 -OUTPUT=0 -GRID=MNGUSG_250m.sgrd -MULTIPLE=0 -LINE_TYPE=1 -POLY_TYPE=1 -GRID_TYPE=1 -TARGET_DEFINITION=0 -TARGET_USER_XMIN=-180 -TARGET_USER_XMAX=180 -TARGET_USER_YMIN=-90 -TARGET_USER_YMAX=90 -TARGET_USER_FITS=1 -TARGET_USER_SIZE=0.002083333'))
#system(paste(gdalwarp, 'MANGPR_250m.sdat  -tr 0.002083333 0.002083333 -ot \"Byte\" -dstnodata 255'))
## TH: must be a bug because it does not expect user defined grid

## manually prepared 500m grid
system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"MANGPR_500m.sgrd\" -FORMULA=\"g1*100\" -INTERPOLATION=0 -USE_NODATA=0 -TYPE=4 -RESULT=\"MANGPRf_500m.sgrd\"'))
unlink("MANGPR_500m.sdat")

## test it:
system(paste(gdalwarp, 'MANGPRf_500m.sdat test.tif -r bilinear -tr 0.002083333 0.002083333 -te -17 10 -14 14 -ot \"Byte\" -dstnodata 255 -overwrite'))
plot(raster("test.tif"), col=SAGA_pal[[1]])
## OK

## Select tiles with Mangroves:
load("/data/models/equi7t1.rda")

## remove all tiles with all missing pixels in the LMK
msk.lst <- list.files(path="/data/covs1t", pattern=glob2rx("MNGUSG_*_*_*.tif"), full.names=TRUE, recursive=TRUE)
length(msk.lst)
i = msk.lst[which(basename(msk.lst)=="MNGUSG_NA_099_020.tif")]

check.LMK <- function(i){
  r <- raster(i)
  na.count <- sum(getValues(is.na(r)|r==0))
  if(na.count==ncell(r)){
    out = 0
  } else {
    if(na.count==0){
      out = 100
    } else {
      out = signif((1-na.count/ncell(r))*100, 3)
    } 
  }
  return(out)
}
sfInit(parallel=TRUE, cpus=48)
sfLibrary(raster)
sfLibrary(rgdal)
sfExport("check.LMK", "msk.lst")
selL <- sfLapply(msk.lst, check.LMK)
sfStop()
## 1692 tiles
summary(unlist(selL))
selD <- data.frame(name=basename(msk.lst)[which(unlist(selL)>0)])
selD$tile <- sapply(paste(selD$name), function(x){strsplit(x, "_")[[1]][2]})
## zip all files:
system(paste('7za a -t7z mangroves_250m.7z ', paste(msk.lst[which(unlist(selL)>0)], collapse=" ")))

## export shapefiles:
for(j in names(equi7t1)){
  sel <- which(equi7t1[[j]]$TILE %in% substr(paste(selD$name[selD$tile==j]), 11, 17))
  out = paste0(j, "_mangroves.shp")
  writeOGR(equi7t1[[j]][sel,], out, ".", "ESRI Shapefile")
}
## compress:
system('7za a -t7z mangroves_shapefile.7z *_mangroves.*')

## Sample data (Key West)
system('7za e N025W085_N030W080.tar.gz')
system('7za e N025W085_N030W080.tar')
ave <- list.files(pattern="AVE_DSM", recursive=TRUE, full.names=TRUE)
cat(ave, sep="\n", file="ave.txt")
system(paste0(gdalbuildvrt, ' -input_file_list ave.txt test.vrt'))

r <- raster(i)
system(paste0(gdalwarp, ' test.vrt DEM_NA_099_020_30m.tif -r bilinear -tr 30 30 -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "),' -t_srs \"', proj4string(r), '\" -overwrite -co \"COMPRESS=DEFLATE\"'))

lnd <- list.files(pattern="B40.TIF.gz", recursive=FALSE, full.names=TRUE)
sapply(lnd, gunzip)
lnd <- list.files(pattern="B40.TIF", recursive=FALSE, full.names=TRUE)
plot(raster(lnd[1]))
cat(lnd, sep="\n", file="lnd.txt")
system(paste0(gdalbuildvrt, ' -input_file_list lnd.txt test2.vrt'))
system(paste0(gdalwarp, ' test2.vrt B40_NA_099_020_30m.tif -r bilinear -tr 30 30 -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "),' -t_srs \"', proj4string(r), '\" -overwrite -co \"COMPRESS=DEFLATE\"'))
