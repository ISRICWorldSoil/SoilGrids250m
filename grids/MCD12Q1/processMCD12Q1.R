## Land cover maps based on MODIS (https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mcd12q1)
## Tom.Hengl@isric.org

library(rgdal)
library(utils)
library(R.utils)
library(snowfall)
library(raster)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")

## Mosaic:
M.lst <- list.files(path="./tiled/", pattern="*.tif$", full.names = TRUE)
N.lst <- list.files(path="./tiledB/", pattern="*.tif$", full.names = TRUE)
GDALinfo(M.lst[1])

make_mosaic <- function(year, typ, lst){
  outm <- paste0("LandCover_", year, "_", typ, "_500m.tif")
  if(!file.exists(outm)){
    lst.s <- lst[grep(pattern=year, lst)]
    vrt <- paste0(typ, '_liste', year, '.vrt')
    txt <- paste0(typ, '_liste', year, '.txt')
    cat(lst.s, sep="\n", file=txt)
    system(paste0(gdalbuildvrt, ' -input_file_list ', txt,' ', vrt))
    system(paste0(gdalwarp, ' ', vrt, ' ', outm, ' -r \"near\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000')) ## -s_srs \"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +units=m +no_defs\" -t_srs \"+proj=longlat +datum=WGS84\"  -tr 0.004166667 0.004166667 -te -180 -90 180 90
    ## ERROR 1: Translating source or target SRS failed
    unlink(vrt)
    unlink(txt)
  }
}
## proper proj4 string: http://gis.stackexchange.com/questions/39116/coordinate-system-of-a-modis-file-to-be-introduced-in-gdal-for-transformation

year.lst = paste0(2001:2013, "001")
sfInit(parallel=TRUE, cpus=48)
sfExport("make_mosaic", "M.lst", "gdalbuildvrt", "gdalwarp", "year.lst")
t <- sfLapply(year.lst, make_mosaic, typ="L1", lst=M.lst)
sfStop()

sfInit(parallel=TRUE, cpus=48)
sfExport("make_mosaic", "N.lst", "gdalbuildvrt", "gdalwarp", "year.lst")
t <- sfLapply(year.lst, make_mosaic, typ="L5", lst=N.lst)
sfStop()

