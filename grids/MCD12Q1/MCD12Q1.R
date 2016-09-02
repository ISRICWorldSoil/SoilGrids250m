## Extract MCD12Q1 land cover band to GeoTiff

library(rgdal)
library(utils)
library(snowfall)
library(raster)
library(RSAGA)
gdal.dir <- shortPathName("C:/Program files/GDAL")
#gdal_setInstallation(search_path=gdal.dir, rescan=TRUE)
#gdal_chooseInstallation(hasDrivers="HDF5")
gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
gdalinfo <- paste0(gdal.dir, "/gdalinfo.exe")
#system(gdal_translate)
modis.prj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
## MODIS Tiles:
load(file="modis_grid.rda")
modis_grid$TILEh <- paste0(ifelse(nchar(modis_grid$h)==1, paste0("h0", modis_grid$h), paste0("h", modis_grid$h)), ifelse(nchar(modis_grid$v)==1, paste0("v0", modis_grid$v), paste0("v", modis_grid$v)))

## Sync all hdf files (RUN AFTER DOWNLOADING WITH FILEZILLA):
y.span = 2001:2013
for(i in y.span){
  dr.lst <- normalizePath(list.dirs(path=paste0("X:\\MODIS\\51\\MCD12Q1\\", i), recursive=FALSE))
  for(j in 1:length(dr.lst)){
    setwd(dr.lst[j])
    x <- strsplit(dr.lst[j], "\\\\")[[1]][6]
    system(paste0('wget --accept \"*.hdf\" -nd -N -r ftp://anonymous@ladsftp.nascom.nasa.gov/allData/51/MCD12Q1/', i ,'/', x)) ## cut-dirs=4
  }
}

hdf.lst <- list.files(path="X:\\MODIS\\51\\MCD12Q1", pattern=glob2rx("*.hdf$"), full.names=TRUE, recursive=TRUE)
## 317 per year (23 folders) --> 354GB 
hdfy.lst <- data.frame(matrix(unlist(lapply(hdf.lst, strsplit, "/")), ncol=4, byrow=T))
hdfy.lst <- hdfy.lst[,-1]
names(hdfy.lst) <- c("Year", "Day", "FileName")
hdfy.lst$YearDay <- as.factor(paste(hdfy.lst$Year, hdfy.lst$Day, sep="_"))
hdfy.lst <- hdfy.lst[hdfy.lst$Year %in% y.span,]
lvs <- levels(as.factor(paste(hdfy.lst$YearDay)))
str(hdfy.lst)
## 4121 tiles
summary(hdfy.lst$YearDay, maxsum=length(levels(hdfy.lst$YearDay)))

setwd("G:/SoilGrids250m/MCD12Q1")
## extract land cover band to GeoTiff:
for(i in 1:length(lvs)){
  sel <- hdfy.lst$YearDay==lvs[i]
  tmp.n <- paste0("./tiled/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  tmp.check <- gsub("./tiled/", "X:/MODIStiled/MCD12Q1/tiled/", tmp.n)
  hddir <- paste0("X:/MODIS/51/MCD12Q1/", hdfy.lst$Year[sel], "/", hdfy.lst$Day[sel], "/")
  tmp0.lst <- paste0('HDF4_EOS:EOS_GRID:\"', hddir, hdfy.lst[sel,"FileName"], '\":MOD12Q1:\"Land_Cover_Type_1\"')
  sfInit(parallel=TRUE, cpus=8)
  sfExport("tmp0.lst", "tmp.check", "tmp.n", "gdal_translate")
  t <- sfLapply(1:length(tmp0.lst), function(j){ if(!file.exists(tmp.check[j])){ try( system(paste(gdal_translate, tmp0.lst[j], tmp.n[j], ' -co \"COMPRESS=DEFLATE\"')) ) }})
  sfStop()
}

## extract also the band "Land Cover Types Description(Type 5 classification)" to GeoTiff:
for(i in 1:length(lvs)){
  sel <- hdfy.lst$YearDay==lvs[i]
  tmp.n <- paste0("./tiledB/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  tmp.check <- gsub("./tiledB/", "X:/MODIStiled/MCD12Q1/tiledB/", tmp.n)
  hddir <- paste0("X:/MODIS/51/MCD12Q1/", hdfy.lst$Year[sel], "/", hdfy.lst$Day[sel], "/")
  tmp0.lst <- paste0('HDF4_EOS:EOS_GRID:\"', hddir, hdfy.lst[sel,"FileName"], '\":MOD12Q1:\"Land_Cover_Type_5\"')
  sfInit(parallel=TRUE, cpus=8)
  sfExport("tmp0.lst", "tmp.check", "tmp.n", "gdal_translate")
  t <- sfLapply(1:length(tmp0.lst), function(j){ if(!file.exists(tmp.check[j])){ try( system(paste(gdal_translate, tmp0.lst[j], tmp.n[j], ' -co \"COMPRESS=DEFLATE\"')) ) }})
  sfStop()
}

