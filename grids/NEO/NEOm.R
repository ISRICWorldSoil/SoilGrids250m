## Average MYDAL2_M_CLD_OT / Cloud Optical Thickness (1 month - Aqua/MODIS)
## Average MYDAL2_M_SKY_WV (WATER VAPOR 1 MONTH - AQUA/MODIS)
## Tom.Hengl@isric.org
## Full description available at: http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MYDAL2_M_CLD_OT
## Effective resolution: 10 km

library(utils)
library(R.utils)
library(snowfall)
library(raster)
gdalwarp =  "/usr/local/bin/gdalwarp"
system("/usr/local/bin/gdal-config --version")
m.lst <- c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
mn.lst <- list(list("01", "02"), list("03", "04"), list("05", "06"), list("07", "08"), list("09", "10"), list("11", "12"))
gdal_translate =  "/usr/local/bin/gdal_translate"

## define functions:
meanf <- function(x){calc(x, mean, na.rm=TRUE)}

for(i in 1:length(mn.lst)){
  M <- c(list.files(path="./MYDAL2_M_SKY_WV", pattern=glob2rx(paste0("*-", mn.lst[[i]][[1]], "-01_rgb_3600x1800.FLOAT.TIFF")), full.names=TRUE), list.files(path="./MYDAL2_M_SKY_WV", pattern=glob2rx(paste0("*-", mn.lst[[i]][[2]], "-01_rgb_3600x1800.FLOAT.TIFF")), full.names=TRUE))
  out.file = paste0("SKY_WV_M_", m.lst[i], "_10km.tif")
  if(!file.exists(out.file)){
    s10km <- raster::stack(M)
    s10km[s10km>100] <- NA 
    ## run in parallel per month:
    beginCluster()
    r <- clusterR(s10km*100, fun=meanf, filename=out.file, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
    endCluster()
  }
}

## filter missing values:
for(i in 1:length(mn.lst)){
  ## check if there are any missing pixels:
  out.file = paste0("SKY_WV_M_", m.lst[i], "_10km.tif")
  r <- raster(out.file)
  if(sum(getValues(is.na(r)))>0){
    system(paste(gdal_translate, out.file, gsub(".tif", ".sdat", out.file), '-ot \"Int16\" -of \"SAGA\"'))
    system(paste0('/usr/local/bin/saga_cmd -c=40 grid_tools 7 -INPUT=\"', gsub(".tif", ".sgrd", out.file), '\" -RESULT=\"tmp.sgrd\"'))
    unlink(out.file)
    system(paste0(gdal_translate, ' tmp.sdat ', out.file, ' -ot \"Int16\" -a_nodata \"-32768\" -co \"COMPRESS=DEFLATE\"'))
    unlink("*.prj")
    unlink("*.sgrd")
    unlink("*.mgrd")
    unlink("*.sdat")
  }
}

save.image("/data/NEO/.RData")
