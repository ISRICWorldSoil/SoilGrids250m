## Extract MOD11A2 LST bands to GeoTiff

library(gdalUtils)
library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
gdal.dir <- shortPathName("C:/Program files/GDAL")
gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe")
#system(gdal_translate)
modis.prj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

y.span = c(2000, 2001, 2004, 2006, 2008, 2011, 2013, 2014)
for(i in y.span){
  dr.lst <- normalizePath(list.dirs(path=paste0("X:\\MODIS\\5\\MOD11A2\\", i), recursive=FALSE))
  for(j in 1:length(dr.lst)){
    setwd(dr.lst[j])
    x <- strsplit(dr.lst[j], "\\\\")[[1]][6]
    system(paste0('wget --accept \"*.hdf\" -nd -N -r ftp://anonymous@ladsftp.nascom.nasa.gov/allData/5/MOD11A2/', i ,'/', x)) ## –cut-dirs=4
  }
}

## list per year
hdf.lst <- list.files(path="X:\\MODIS\\5\\MOD11A2", pattern=glob2rx("*.hdf$"), full.names=TRUE, recursive=TRUE)
hdfy.lst <- data.frame(matrix(unlist(lapply(hdf.lst, strsplit, "/")), ncol=4, byrow=T))
hdfy.lst <- hdfy.lst[,-1]
names(hdfy.lst) <- c("Year", "Day", "FileName")
hdfy.lst$YearDay <- as.factor(paste(hdfy.lst$Year, hdfy.lst$Day, sep="_"))
#x <- gdalinfo(hdf.lst[i])
hdfy.lst <- hdfy.lst[hdfy.lst$Year %in% y.span,]
lvs <- levels(as.factor(paste(hdfy.lst$YearDay)))
str(hdfy.lst)
## 98,616 tiles

## Daytime LST:
for(i in 1:length(lvs)){
  sel <- hdfy.lst$YearDay==lvs[i]
  tmp.n <- paste0("X:/MODIStiled/MOD11A2/tiledD/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  out.n <- paste0("G:/SoilGrids250m/MOD11A2/tiledD/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  #outn <- paste0("MOD11A2_LSTD_", lvs[i], "_1km")
  hddir <- paste0("X:/MODIS/5/MOD11A2/", hdfy.lst$Year[sel], "/", hdfy.lst$Day[sel], "/")
  tmp0.lst <- paste0('HDF4_EOS:EOS_GRID:\"', hddir, hdfy.lst[sel,"FileName"], '\":MODIS_Grid_8Day_1km_LST:\"LST_Day_1km\"') 
  sfInit(parallel=TRUE, cpus=8)
  sfExport("tmp0.lst", "tmp.n", "out.n", "gdal_translate")
  t <- sfLapply(1:length(tmp0.lst), function(j){ if(!file.exists(tmp.n[j])){ try( system(paste(gdal_translate, tmp0.lst[j], out.n[j], ' -a_nodata \"0\" -ot \"Int16\" -co \"COMPRESS=DEFLATE\"')) ) }})
  sfStop()
}

## Nightime LST:
for(i in 1:length(lvs)){
  sel <- hdfy.lst$YearDay==lvs[i]
  tmp.n <- paste0("X:/MODIStiled/MOD11A2/tiledN/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  out.n <- paste0("G:/SoilGrids250m/MOD11A2/tiledN/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  #outn <- paste0("MOD11A2_LSTN_", lvs[i], "_1km")
  hddir <- paste0("X:/MODIS/5/MOD11A2/", hdfy.lst$Year[sel], "/", hdfy.lst$Day[sel], "/")
  tmp0.lst <- paste0('HDF4_EOS:EOS_GRID:\"', hddir, hdfy.lst[sel,"FileName"], '\":MODIS_Grid_8Day_1km_LST:\"LST_Night_1km\"') 
  sfInit(parallel=TRUE, cpus=8)
  sfExport("tmp0.lst", "tmp.n", "out.n", "gdal_translate")
  t <- sfLapply(1:length(tmp0.lst), function(j){ if(!file.exists(tmp.n[j])){ try( system(paste(gdal_translate, tmp0.lst[j], out.n[j], ' -a_nodata \"0\" -ot \"Int16\" -co \"COMPRESS=DEFLATE\"')) ) }})
  sfStop()
}

## Output is 2 x 70GB
#system('7za a -r MOD11A2tiles.zip G:\\SoilGrids250m\\MOD11A2', show.output.on.console=FALSE)