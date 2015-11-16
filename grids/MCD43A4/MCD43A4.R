## Extract MCD43A4 band 4 and 7 to GeoTiff

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
load("modis_grid.rda")

## Sync all hdf files (AFTER DOWNLOADING WITH FILEZILLA):
for(i in c(2006, 2010, 2014)){
  dr.lst <- normalizePath(list.dirs(path=paste0("X:\\MODIS\\5\\MCD43A4\\", i), recursive=FALSE))
  for(j in 1:length(dr.lst)){
    setwd(dr.lst[j])
    x <- strsplit(dr.lst[j], "\\\\")[[1]][6]
    system(paste0('wget --accept \"*.hdf\" -nd -N -r ftp://anonymous@ladsftp.nascom.nasa.gov/allData/5/MCD43A4/', i ,'/', x)) ## –cut-dirs=4
  }
}

## 14,691 file per year = 39,973 files (23 folders per year) --> 850GB 
## list per year
hdf.lst <- list.files(path="X:\\MODIS\\5\\MCD43A4", pattern=glob2rx("*.hdf$"), full.names=TRUE, recursive=TRUE)
hdfy.lst <- data.frame(matrix(unlist(lapply(hdf.lst, strsplit, "/")), ncol=4, byrow=T)
)
hdfy.lst <- hdfy.lst[,-1]
names(hdfy.lst) <- c("Year", "Day", "FileName")
hdfy.lst$YearDay <- as.factor(paste(hdfy.lst$Year, hdfy.lst$Day, sep="_"))
lvs <- levels(hdfy.lst$YearDay)
#x <- gdalinfo(hdf.lst[i])
str(hdfy.lst)
## 18,500 tiles

## extract NRB 4 to GeoTiff:
for(i in 1:length(lvs)){
  sel <- hdfy.lst$YearDay==lvs[i]
  tmp.n <- paste0("./tiled4/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  hddir <- paste0("X:/MODIS/5/MCD43A4/", hdfy.lst$Year[sel], "/", hdfy.lst$Day[sel], "/")
  tmp0.lst <- paste0('HDF4_EOS:EOS_GRID:\"', hddir, hdfy.lst[sel,"FileName"], '\":MOD_Grid_BRDF:\"Nadir_Reflectance_Band4\"')
  sfInit(parallel=TRUE, cpus=8)
  sfExport("tmp0.lst", "tmp.n", "gdal_translate")
  t <- sfLapply(1:length(tmp0.lst), function(j){ if(!file.exists(tmp.n[j])){ try( system(paste(gdal_translate, tmp0.lst[j], tmp.n[j], '-co \"COMPRESS=DEFLATE\"')) ) }})
  sfStop()
}
#outn <- paste0("MCD43A4_NRB4_", lvs[i], "_500m")

## extract NRB 7 to GeoTiff:
for(i in 1:length(lvs)){
  sel <- hdfy.lst$YearDay==lvs[i]
  tmp.n <- paste0("./tiled7/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  hddir <- paste0("X:/MODIS/5/MCD43A4/", hdfy.lst$Year[sel], "/", hdfy.lst$Day[sel], "/")
  tmp0.lst <- paste0('HDF4_EOS:EOS_GRID:\"', hddir, hdfy.lst[sel,"FileName"], '\":MOD_Grid_BRDF:\"Nadir_Reflectance_Band7\"')
  sfInit(parallel=TRUE, cpus=8)
  sfExport("tmp0.lst", "tmp.n", "gdal_translate")
  t <- sfLapply(1:length(tmp0.lst), function(j){ if(!file.exists(tmp.n[j])){ try( system(paste(gdal_translate, tmp0.lst[j], tmp.n[j], '-co \"COMPRESS=DEFLATE\"')) ) }})
  sfStop()
}

#outn <- paste0("MCD43A4_NRB7_", lvs[i], "_500m")
load("dirs.rda")
Tiles <- expand.grid(modis_grid$TILEh, dirs)
et1 <- sapply(list.files(path="X:/MODIS/5/MCD43A4/2014/025", pattern=glob2rx("*.hdf$")), function(x){strsplit(x, "\\.")[[1]][3]})
et2 <- sapply(list.files(path="X:/MODIS/5/MCD43A4/2006/345", pattern=glob2rx("*.hdf$")), function(x){strsplit(x, "\\.")[[1]][3]}) 
et3 <- sapply(list.files(path="X:/MODIS/5/MCD43A4/2006/345", pattern=glob2rx("*.hdf$")), function(x){strsplit(x, "\\.")[[1]][3]})
tiles.n <- levels(as.factor(paste(c(et1,et2,et3))))
## 322 tiles

modis_grid_MCD43A4.xy <- data.frame(TILEh=tiles.n, ulx=rep(NA, length(tiles.n)), uly=rep(NA, length(tiles.n)), lrx=rep(NA, length(tiles.n)), lry=rep(NA, length(tiles.n)))
for(k in 1:nrow(modis_grid_MCD43A4.xy)){
  try( r <- raster(list.files(path="G:\\SoilGrids250m\\MCD43A4\\tiled4", pattern=glob2rx(paste0("MCD43A4.A2014025.", tiles.n[k],".*.*.tif$")), full.names=TRUE)) )
  try( modis_grid_MCD43A4.xy[k,c("ulx","uly","lrx","lry")] <- as.vector(extent(r))[c(1,4,2,3)] )
}
save(modis_grid_MCD43A4.xy, file="modis_grid_MCD43A4.xy.rda")
    
    