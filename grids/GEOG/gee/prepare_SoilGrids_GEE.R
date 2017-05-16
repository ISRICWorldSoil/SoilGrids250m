## Prepare SoilGrids for import into GEE
## by tom.hengl@isric.org and milan.kili11@gmail.com

library(rgdal)
library(snowfall)
var.tbl = read.csv("list_gge_names.csv")
str(var.tbl)

bind_rename <- function(var, gee.name, bands){
  in.files = list.files(path="/data/GEOG", pattern=glob2rx(paste0(var, "*.tif$")), full.names = TRUE)
  if(length(in.files)>1 & bands>1){
    if(!file.exists(paste0(gee.name, '.tif'))){
      out.tmp <- tempfile(fileext = ".txt")
      cat(in.files, sep="\n", file=out.tmp)
      vrt = paste0(var, ".vrt")
      system(paste0('gdalbuildvrt -separate ', vrt, ' -input_file_list ', out.tmp))
      system(paste0('gdal_translate --config GDAL_CACHEMAX 99000 ', vrt,' ', gee.name, '.tif -co \"COMPRESS=DEFLATE\" -co BIGTIFF=YES'))
    }
  } else {
    system(paste0('cp ', in.files, ' /data/GEOG/gee/', gee.name, '.tif'))
  }
}

sfInit(parallel=TRUE, cpus=5)
sfExport("var.tbl", "bind_rename")
sfLibrary(rgdal)
out <- sfClusterApplyLB(1:nrow(var.tbl), function(i){bind_rename(var=paste(var.tbl[i,"SoilGrids_code"]), gee.name=paste(var.tbl[i,"Asset_name"]), bands=var.tbl[i,"Bands"])})
sfStop()
