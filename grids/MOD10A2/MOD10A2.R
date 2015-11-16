## Extract MOD10A2 to GeoTiff

library(gdalUtils)
library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
gdal.dir <- shortPathName("C:/Program files/GDAL")
gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)
#gdal_chooseInstallation(hasDrivers="HDF5")
gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
#system(gdal_translate)
modis.prj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

## Sync all hdf files (AFTER DOWNLOADING WITH FILEZILLA):
dr.lst <- normalizePath(list.dirs(path="X:\\MODIS\\5\\MOD10A2", recursive=FALSE))
for(j in 1:length(dr.lst)){
  setwd(dr.lst[j])
  x <- strsplit(dr.lst[j], "\\\\")[[1]][5]
  system(paste0('wget --accept \"*.hdf\" -nd -N -r ftp://n5eil01u.ecs.nsidc.org/SAN/MOST/MOD10A2.005/', x))
}

## 152,048 hdf files --> 35GB 
## TAKES CA 10 MINS TO LIST ALL FILES!!
hdf.lst <- list.files(path="X:\\MODIS\\5\\MOD10A2", pattern=glob2rx("*.hdf$"), full.names=TRUE, recursive=TRUE)
hdfy.lst <- data.frame(matrix(unlist(lapply(hdf.lst, strsplit, "/")), ncol=3, byrow=T))
hdfy.lst <- hdfy.lst[,-1]
names(hdfy.lst) <- c("Year_Day", "FileName")
lvs <- levels(as.factor(hdfy.lst$Year_Day))
## 345
#x <- gdalinfo(hdf.lst[i])
str(hdfy.lst)
## focus only on complete years:
lvs <- lvs[!sapply(lvs, function(x){strsplit(x, "\\.")[[1]][1]})==2015]
## 321

setwd("G:/SoilGrids250m/MOD10A2")
## extract Eight_Day_Snow_Cover to GeoTiff:
for(i in 1:length(lvs)){
  sel <- hdfy.lst$Year_Day==lvs[i]
  tmp.n <- paste0("./tiled/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  hddir <- paste0("X:/MODIS/5/MOD10A2/", hdfy.lst$Year[sel], "/", hdfy.lst$Day[sel], "/")
  tmp0.lst <- paste0('HDF4_EOS:EOS_GRID:\"', hddir, hdfy.lst[sel,"FileName"], '\":MOD_Grid_Snow_500m:\"Eight_Day_Snow_Cover\"')
  sfInit(parallel=TRUE, cpus=8)
  sfExport("tmp0.lst", "tmp.n", "gdal_translate")
  t <- sfLapply(1:length(tmp0.lst), function(j){ if(!file.exists(tmp.n[j])){ try( system(paste(gdal_translate, tmp0.lst[j], tmp.n[j], '-co \"COMPRESS=DEFLATE\"')) ) }})
  sfStop()
}
## 97,720 tiles for years 2001, 2002, 2010, 2011, 2012, 2013, 2014!

#outn <- paste0("MOD10A2_Snow_", lvs[i], "_500m")

    
    