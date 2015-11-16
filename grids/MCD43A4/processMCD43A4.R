## Average RB4 and RB7 at 500 m
## by: Tom.Hengl@isric.org
## Example tile name:
## MCD43A4.A2010081.h08v05.005.2010101144907.hdf
## each time is 4800 x 4800 pixels
## 20,196 input tiles for 3 years

library(utils)
library(R.utils)
library(RSAGA)
library(snowfall)
library(raster)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
system('/usr/local/bin/saga_cmd --version')

## MODIS tiles and days of the year
load(file="modis_grid.rda")
str(modis_grid$TILE)
## only 317 cover land
modis_grid$TILEh <- paste0(ifelse(nchar(modis_grid$h)==1, paste0("h0", modis_grid$h), paste0("h", modis_grid$h)), ifelse(nchar(modis_grid$v)==1, paste0("v0", modis_grid$v), paste0("v", modis_grid$v)))
## 8-day period
load("dirs.rda")
days <- as.numeric(format(seq(ISOdate(2015,1,1), ISOdate(2015,12,31), by="month"), "%j"))-1
m.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dsel <- cut(as.numeric(dirs), breaks=c(days,365), labels = m.lst)
modis.prj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
xs <- list.files(path="./tilednir", pattern=glob2rx("*.A2014257.*.*.*.tif$"))
## 326 x 12 MODIS tiles

aggr_tile <- function(lst, outm){
  if(!file.exists(outm)){
    s <- raster::stack(lst)
    smean <- raster::calc(s, fun=mean, na.rm=TRUE)*10000
    r1 <- raster::writeRaster(smean, filename=outm, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
  }
}

aggr_parallel <- function(i, path, channel, band){
  lstD <- list.files(path=paste0("./tiled", channel), pattern=glob2rx(paste0("MCD43A4.*.", modis_grid$TILEh[i], ".005.*.tif$")), full.names=TRUE)
  if(length(lstD)>0){
    for(k in 1:length(m.lst)){
      lstD.s <- lstD[unlist(sapply(dirs[which(m.lst[k] == dsel)], function(x){grep(basename(lstD), pattern=glob2rx(paste0("*.A2???",x,".*.*.*.tif$")))}))]
      if(length(lstD.s)>0){
        try( aggr_tile(lstD.s, outm=paste0("./Mtiled", channel, "/NRB", band, "_M_", m.lst[k], "_", modis_grid$TILEh[i], ".tif")) )
      }
    }
  }
}

## test it:
aggr_parallel(i=421, channel="nir", band=4)
aggr_parallel(i=565, channel="mir", band=7)

sfInit(parallel=TRUE, cpus=35)
sfExport("modis_grid", "aggr_tile", "m.lst", "dsel", "dirs")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(sp)
t <- sfLapply(1:length(modis_grid$TILEh), aggr_parallel, channel="nir", band=4)
sfStop()

sfInit(parallel=TRUE, cpus=35)
sfExport("modis_grid", "aggr_tile", "m.lst", "dsel", "dirs")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(sp)
t <- sfLapply(1:length(modis_grid$TILEh), aggr_parallel, channel="mir", band=7)
sfStop()

## Mosaic:
M4.lstD <- list.files(path="./Mtilednir", pattern=glob2rx("NRB4_M_*_*.tif$"), full.names = TRUE)
#v01.rm <- M4.lstD[grep(pattern="v01", M4.lstD)]
M7.lstD <- list.files(path="./Mtiledmir", pattern=glob2rx("NRB7_M_*_*.tif$"), full.names = TRUE)
## 3719 tiles x 2

make_mosaic <- function(j, band){
  outm <- paste0("NRB", band, "_M_", m.lst[j], "_500m.tif")
  if(!file.exists(paste0(outm,'.gz'))){
    lst <- get(paste0("M", band, ".lstD"))
    lst.s <- lst[grep(pattern=m.lst[j], lst)]
    vrt <- paste0('M_liste', m.lst[j], '_', band, '.vrt')
    txt <- paste0('M_liste', m.lst[j], '_', band, '.txt')
    cat(lst.s, sep="\n", file=txt)
    system(paste0(gdalbuildvrt, ' -input_file_list ', txt,' ', vrt))
    system(paste0(gdalwarp, ' ', vrt, ' ', outm, ' -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Int16\" -dstnodata \"-32768\" -tr 0.004166667 0.004166667 -te -180 -90 180 90'))
    gzip(outm)
    #unlink(vrt)
    #unlink(txt)
  }
}

## test:
make_mosaic(j=2, band=7)

bands <- c(rep(4,12), rep(7,12))
mon <- rep(1:length(m.lst), 2)

sfInit(parallel=TRUE, cpus=12)
sfExport("make_mosaic", "M4.lstD", "M7.lstD", "m.lst", "gdalbuildvrt", "gdalwarp", "mon", "bands")
sfLibrary(R.utils)
sfLibrary(sp)
t <- sfLapply(1:length(mon), function(x){ make_mosaic(mon[x], band=bands[x])})
sfStop()

## TH: Some months miss many pixels, especially Jan and Feb for latitudes above 63 degree.
## TH: soilution is to smooth out values using the neighbouring months
## Mean annual temperature:
Ml <- paste0("NRB4_M_", m.lst, "_500m.tif")
sapply(paste0(Ml,".gz"), gunzip)
s500m <- raster::stack(Ml)
## define functions:
smoothf <- function(x){calc(x, mean, na.rm=TRUE)}
## run in parallel:
beginCluster(n=40)
for(i in 1:length(m.lst)){
    if(i==1){ s <- raster::stack(raster(s500m, 11), raster(s500m, 12), raster(s500m, 1), raster(s500m, 1), raster(s500m, 1), raster(s500m, 2), raster(s500m, 3))  }
    if(i==12){ s <- raster::stack(raster(s500m, 10), raster(s500m, 11), raster(s500m, 12), raster(s500m, 12), raster(s500m, 12), raster(s500m, 1), raster(s500m, 2))  }
    if(i>1&i<12){ s <- raster::stack(raster(s500m, i-1), raster(s500m, i), raster(s500m, i), raster(s500m, i+1)) }
    outn = paste0("NRB4_Ms_", m.lst[i], "_500m.tif")
    if(!file.exists(outn)){  
      r <- clusterR(s, fun=smoothf, filename=outn, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
    }
}
endCluster()

