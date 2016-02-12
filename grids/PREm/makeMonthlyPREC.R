## Derivation of mean monthly precipitation using 2 data sources:
## WorldClim (http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/prec_30s_bil.zip) and GPCP Combined Precipitation Dataset (ftp://ftp.cdc.noaa.gov/Datasets/gpcc/full_v6/precip.mon.1981-2010.ltm.v6.nc)
## NASA's GPCP project is at http://precip.gsfc.nasa.gov/

library(sp)
library(rgdal)
library(raster)
gdalwarp = "/usr/local/bin/gdalwarp"
gdal_translate = "/usr/local/bin/gdal_translate"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
gdaladdo = "/usr/local/bin/gdaladdo"
system("/usr/local/bin/gdal-config --version")
#GDALinfo("prec_1_10km.sdat")
#plot(raster("prec_1_10km.sdat"))

unzip("GPCC_monthly.zip") ## GPCP LTM Total Full V6 (0.5x0.5 degree)
unzip("prec_30s_bil.zip") ## WorldClim 1km
bil.lst <- list.files(pattern=".bil$")
## Convert to SAGA GIS format:
trans <- function(x){
  out = gsub(".bil", ".sdat", x)
  if(!file.exists(out)){
    system(paste0(gdalwarp, ' ', x, ' ',out, ' -of \"SAGA\" -te -180 -90 180 90 -tr 0.008333333 0.008333333'))
  }
}
sfInit(parallel=TRUE, cpus=12)
sfExport("trans", "bil.lst", "gdalwarp")
sfLibrary(rgdal)
x <- sfClusterApplyLB(bil.lst, trans)
sfStop()
GDALinfo("prec_1.sdat")

lst.gpcp <- paste0("prec_", 1:12, "_10km.sgrd")
lst.wc <- paste0("prec_", 1:12, ".sgrd")
mon <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
cod <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")


## Derive average of two maps - WorldClim and GPCP (satellite based PREC)
## define functions:
meanf <- function(x){calc(x, mean, na.rm=TRUE)}
sumf <- function(x){calc(x, sum, na.rm=TRUE)}

## TAKES CA 1hr
for(j in 1:length(lst.gpcp)){
  outtif <- paste0("P", cod[j], "MRG3a.tif")
  if(!file.exists(outtif)){
    if(!file.exists(paste0("PREm_", j, "_1km.sgrd"))){
      ## downscale to 1 km:
      system(paste0('/usr/local/bin/saga_cmd -c=48 grid_tools 0 -INPUT=\"', lst.gpcp[j], '\" -SCALE_DOWN=3 -TARGET_DEFINITION=0 -TARGET_USER_XMIN=', -180+1/240,' -TARGET_USER_XMAX=', 180-1/240,' -TARGET_USER_YMIN=', -60+1/240,' -TARGET_USER_YMAX=', 90-1/240,' -TARGET_USER_SIZE=0.0083333333 -OUTPUT=\"PREm_', j, '_1km.sgrd\"'))
      #GDALinfo("PREm_1_1km.sdat")
    }
    ## calculate average:
    sN <- raster::stack(c(gsub(".sgrd", ".sdat", lst.wc[j]), paste0('PREm_', j, '_1km.sdat'), paste0('PREm_', j, '_1km.sdat'))) ## two times more weight to GPCP!
    ## run in parallel:
    beginCluster()
    r <- clusterR(sN, fun=meanf, filename=outtif, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
    endCluster()
    system(paste0(gdaladdo, ' ', outtif, ' 2 4 8 16 32 64'))
  }
}

## clean-up:
unlink(list.files(pattern=".sdat$"))
unlink(list.files(pattern=".sgrd$"))
unlink(list.files(pattern=".mgrd$"))
unlink(list.files(pattern=".prj$"))
unlink(list.files(pattern=".bil$"))
unlink(list.files(pattern=".hdr$"))

## Mean annual precipitation:
sN <- raster::stack(list.files(pattern="MRG3a"))
beginCluster()
r <- clusterR(sN, fun=sumf, filename="PRSMRG3a.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
endCluster()
system(paste0(gdaladdo, ' PRSMRG3a.tif 2 4 8 16 32 64'))

