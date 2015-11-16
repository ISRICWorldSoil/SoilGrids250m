## Long term Average EVI at 250 m
## Tom.Hengl@isric.org
## Example tile name:
## MOD13Q1.A2006001.h12v04.005.2008063151403.tif
## each tile is 4800 x 4800 pixels
## 40,392 input tiles for 7 years

library(utils)
library(R.utils)
library(RSAGA)
library(snowfall)
library(raster)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")

## MODIS tiles and days of the year
load(file="modis_grid.rda")
str(modis_grid$TILE)
## only 317 cover land
modis_grid$TILEh <- paste0(ifelse(nchar(modis_grid$h)==1, paste0("h0", modis_grid$h), paste0("h", modis_grid$h)), ifelse(nchar(modis_grid$v)==1, paste0("v0", modis_grid$v), paste0("v", modis_grid$v)))
## 8-day period
load("dirs.rda")
days <- as.numeric(format(seq(ISOdate(2015,1,1), ISOdate(2015,12,31), by="month"), "%j"))-1
#m.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
## we use 6 periods to get smoother values / less missing values
## soils are not per correlated with current biomass / land cover
m.lst <- c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
dsel <- cut(as.numeric(dirs), breaks=c(days[c(1,3,5,7,9,11)], 365), labels = m.lst)
modis.prj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m+no_defs"

## Parallel processing with SAGA GIS:
aggr_tile <- function(lst, outm, outs){
  if(!file.exists(outm)){
    ## convert to SAGA format:
    tmp <- list(NULL)
    for(j in 1:length(lst)){
      tmp[[j]] <- set.file.extension(tempfile(tmpdir="./Mtiled"), ".sdat")
      try( system(paste0(gdal_translate, ' ', lst[j], ' ', tmp[[j]], ' -q -of \"SAGA\" -ot \"Int16\" -tr 231.656 231.656 -r \"near\"')) )
    }
    ## derive average and STD using all 40 cores:
    try( system(paste0('/usr/local/bin/saga_cmd -c=38 statistics_grid 4 -GRIDS \"', paste(set.file.extension(unlist(tmp), ".sgrd"), collapse=";", sep=""), '\" -MEAN ', outm, ' -STDDEV ', outs)) )
    unlink(tmp)
    unlink(set.file.extension(unlist(tmp), ".prj"))
    unlink(set.file.extension(unlist(tmp), ".sgrd"))
    unlink(set.file.extension(unlist(tmp), ".mgrd"))
  }  
}

## Generate average per month (per year):
aggr_parallel <- function(i){
  lstD <- list.files(path="./tiled/", pattern=glob2rx(paste0("MOD13Q1.*.", modis_grid$TILEh[i], ".005.*.tif$")), full.names=TRUE)
  if(length(lstD)>0){
    for(k in 1:length(m.lst)){
      lstD.s <- lstD[unlist(sapply(dirs[which(m.lst[k] == dsel)], function(x){grep(basename(lstD), pattern=glob2rx(paste0("*.A2???",x,".*.*.*.tif$")))}))] 
      if(length(lstD.s)>0){
        aggr_tile(lstD.s, outm=paste0("./Mtiled/EVI_M_", m.lst[k], "_", modis_grid$TILEh[i], ".sgrd"), outs=paste0("./Mtiled/EVI_SD_", m.lst[k], "_", modis_grid$TILEh[i], ".sgrd"))
      }
    }
  }
}

t <- sapply(1:length(modis_grid$TILEh), aggr_parallel)
## TAKES ca 19 HRS!
## output -> 296 tiles * 6 periods * 2 = 3552 images
## cleanup
unlink(list.files(path="./Mtiled", pattern="file", full.names=TRUE))

## Convert back to GeoTiff and fix coordinates:
sdat.lst <- list.files(path="./Mtiled/", pattern="*.sdat$", full.names = TRUE)

sdat2gdal <- function(j, sdat.lst){
  try( system(paste0(gdal_translate, ' ', sdat.lst[j], ' ', set.file.extension(sdat.lst[j], ".tif"), ' -ot \"Int16\" -a_nodata \"-32768\" -a_srs \"', modis.prj,'\" -co \"COMPRESS=DEFLATE\"')) )
  if(file.exists(set.file.extension(sdat.lst[j], ".tif"))){
    unlink(sdat.lst[j])
    unlink(set.file.extension(sdat.lst[j], ".prj"))
    unlink(set.file.extension(sdat.lst[j], ".sgrd"))
    unlink(set.file.extension(sdat.lst[j], ".mgrd"))
  }
}

sfInit(parallel=TRUE, cpus=38)
sfExport("sdat.lst", "sdat2gdal", "gdal_translate", "modis.prj")
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(RSAGA)
t4 <- sfLapply(1:length(sdat.lst), sdat2gdal, sdat.lst)
sfStop()

## Mosaic:
M.lstD <- list.files(path="./Mtiled/", pattern=glob2rx("EVI_M_*_*.tif$"), full.names = TRUE)
SD.lstD <- list.files(path="./Mtiled/", pattern=glob2rx("EVI_SD_*_*.tif$"), full.names = TRUE)
## 3552 tiles

make_mosaic <- function(j, mon, typ){
  outm <- paste0("EVI", "_", typ[j], "_", mon[j], "_250m.tif")
  if(!file.exists(paste0(outm,'.gz'))){
    lst <- get(paste0(typ[j],".lstD"))
    lst.s <- lst[grep(pattern=mon[j], lst)]
    vrt <- paste0(typ[j], '_liste', mon[j], '.vrt')
    txt <- paste0(typ[j], '_liste', mon[j], '.txt')
    cat(lst.s, sep="\n", file=txt)
    system(paste0(gdalbuildvrt, ' -input_file_list ', txt,' ', vrt))
    system(paste0(gdalwarp, ' ', vrt, ' ', outm, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"bilinear\" -ot \"Int16\" -dstnodata \"-32768\" -tr 0.002083333 0.002083333 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))
    #gzip(outm)
    #unlink(vrt)
    #unlink(txt)
  }
}

mon <- rep(m.lst, 2)
typ <- as.vector(sapply(c("M","SD"), rep, length(m.lst)))

sfInit(parallel=TRUE, cpus=12)
sfExport("make_mosaic", "M.lstD", "SD.lstD", "m.lst", "gdalbuildvrt", "gdalwarp", "mon", "typ")
sfLibrary(R.utils)
sfLibrary(sp)
t <- sfLapply(1:length(mon), make_mosaic, mon=mon, typ=typ)
sfStop()
