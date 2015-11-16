## Extract MOD13Q1 EVI band to GeoTiff

library(gdalUtils)
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
#system(gdal_translate)
modis.prj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
## MODIS Tiles:
load(file="modis_grid.rda")
modis_grid$TILEh <- paste0(ifelse(nchar(modis_grid$h)==1, paste0("h0", modis_grid$h), paste0("h", modis_grid$h)), ifelse(nchar(modis_grid$v)==1, paste0("v0", modis_grid$v), paste0("v", modis_grid$v)))
load("dirs.rda")
Tiles <- expand.grid(modis_grid$TILEh, dirs)

## Sync all hdf files (AFTER DOWNLOADING WITH FILEZILLA):
y.span = c(2006, 2008, 2009, 2010, 2012, 2014)
for(i in y.span){
  dr.lst <- normalizePath(list.dirs(path=paste0("X:\\MODIS\\5\\MOD13Q1\\", i), recursive=FALSE))
  for(j in 1:length(dr.lst)){
    setwd(dr.lst[j])
    x <- strsplit(dr.lst[j], "\\\\")[[1]][6]
    system(paste0('wget --accept \"*.hdf\" -nd -N -r ftp://anonymous@ladsftp.nascom.nasa.gov/allData/5/MOD13Q1/', i ,'/', x)) ## –cut-dirs=4
  }
}

## ?? files per year (23 folders) --> 610GB 
## list per year
hdf.lst <- list.files(path="X:\\MODIS\\5\\MOD13Q1", pattern=glob2rx("*.hdf$"), full.names=TRUE, recursive=TRUE)
hdfy.lst <- data.frame(matrix(unlist(lapply(hdf.lst, strsplit, "/")), ncol=4, byrow=T))
hdfy.lst <- hdfy.lst[,-1]
names(hdfy.lst) <- c("Year", "Day", "FileName")
hdfy.lst$YearDay <- as.factor(paste(hdfy.lst$Year, hdfy.lst$Day, sep="_"))
hdfy.lst <- hdfy.lst[hdfy.lst$Year %in% y.span,]
lvs <- levels(as.factor(paste(hdfy.lst$YearDay)))
str(hdfy.lst)
## 40,392 tiles
summary(hdfy.lst$YearDay, maxsum=length(levels(hdfy.lst$YearDay)))

## extract EVI band to GeoTiff:
for(i in 1:length(lvs)){
  sel <- hdfy.lst$YearDay==lvs[i]
  tmp.n <- paste0("./tiled/", sapply(paste(hdfy.lst[sel,"FileName"]), function(x){strsplit(x, ".hdf")[[1]][1]}), ".tif")
  tmp.check <- gsub("./tiled/", "X:/MODIStiled/MOD13Q1/tiled/", tmp.n)
  outn <- paste0("MOD13Q1_EVI_", lvs[i], "_250m")
  hddir <- paste0("X:/MODIS/5/MOD13Q1/", hdfy.lst$Year[sel], "/", hdfy.lst$Day[sel], "/")
  tmp0.lst <- paste0('HDF4_EOS:EOS_GRID:\"', hddir, hdfy.lst[sel,"FileName"], '\":MODIS_Grid_16DAY_250m_500m_VI:\"250m 16 days EVI\"')
  sfInit(parallel=TRUE, cpus=8)
  sfExport("tmp0.lst", "tmp.check", "tmp.n", "gdal_translate")
  t <- sfLapply(1:length(tmp0.lst), function(j){ if(!file.exists(tmp.check[j])){ try( system(paste(gdal_translate, tmp0.lst[j], tmp.n[j], '-ot \"Int16\" -co \"COMPRESS=DEFLATE\"')) ) }})
  sfStop()
}

## Check if all stack together:
Tiles[Tiles$Var1=="h22v13"&Tiles$Var2=="017",]

check_extent <- function(i){
  nm <- list.files(path="./tiled", pattern=glob2rx(paste0("*.A2???",Tiles$Var2[i],".", Tiles$Var1[i],".*.*.tif$")), full.names=TRUE)
  if(length(nm)>0){
    gd <- lapply(nm, GDALinfo, silent=TRUE)
    ll.x <- sapply(gd, function(x){x[4]})
    ll.y <- sapply(gd, function(x){x[5]})
    if(length(unique(ll.x))==1&length(unique(ll.y))==1){
      out <- "OK"
    } else {
      out <- paste0("A2???",Tiles$Var2[i],".",Tiles$Var1[i])
    }
  } else {
    out <- "missing"
  } 
  return(out)
}

## TAKES CA 15 MINS TO CHECK
sfInit(parallel=TRUE, cpus=10)
sfExport("Tiles", "check_extent")
sfLibrary(sp)
sfLibrary(rgdal)
Tiles$check <- unlist(sfLapply(1:nrow(Tiles), check_extent))
sfStop()
summary(as.factor(Tiles$check))

Tiles.complete <- Tiles[!Tiles$check=="missing",]
write.csv(Tiles.complete, "Tiles.complete.csv")
str(Tiles.complete)
## 6732 levels

sel.del <- Tiles$check[!(Tiles$check=="missing"|Tiles$check=="OK")]
for(k in 1:length(sel.del)){
  unlink(list.files(path="./tiled", pattern=glob2rx(paste0("*.", sel.del[k],".*.*.tif$")), full.names=TRUE))
}

## Output is 244GB!!
## CAN TAKE 1 days just to copy to server
#system('7za a -r MOD13Q1tiled.zip G:\\SoilGrids250m\\MOD13Q1\\tiled', show.output.on.console=FALSE)

tiles.n <- levels(as.factor(paste(Tiles.complete$Var1)))
## 296 tiles
modis_grid.xy <- data.frame(TILEh=tiles.n, ulx=rep(NA, length(tiles.n)), uly=rep(NA, length(tiles.n)), lrx=rep(NA, length(tiles.n)), lry=rep(NA, length(tiles.n)))
for(k in 1:nrow(modis_grid.xy)){
  r <- raster(list.files(path="./tiled", pattern=glob2rx(paste0("*.A2010129.", tiles.n[k],".*.*.tif$")), full.names=TRUE))
  modis_grid.xy[k,c("ulx","uly","lrx","lry")] <- as.vector(extent(r))[c(1,4,2,3)]
}
save(modis_grid.xy, file="modis_grid.xy.rda")

#sdat.lst <- list.files(path="G:/BEAR/MOD13Q1/Mtiled/", pattern="*.sdat$", full.names = TRUE)
#sdat2gdal <- function(k){
#  tn <- strsplit(strsplit(sdat.lst[k], "_")[[1]][4], ".sdat")[[1]][1]
#  try( system(paste0(gdal_translate, ' ', sdat.lst[k], ' ', set.file.extension(sdat.lst[k], ".tif"), ' -ot \"Int16\" -a_ullr ', paste(modis_grid.xy[modis_grid.xy$TILEh==tn,c("ulx","uly","lrx","lry")], collapse=" "), ' -a_nodata \"-32768\" -a_srs \"', modis.prj,'\"')) ) ## -tr 231.65635826375006 231.65635826375006
#}    
    