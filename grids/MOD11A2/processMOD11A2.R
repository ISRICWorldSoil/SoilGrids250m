## Average LST Day and Night time
## Tom.Hengl@isric.org
## Example tile name:
## MOD11A2.A2001001.h02v06.005.2006352014912.tif

library(utils)
library(R.utils)
library(snowfall)
library(raster)
gdal_translate =  "/usr/local/bin/gdal_translate"
gdalwarp =  "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
system("/usr/local/bin/gdal-config --version")
## remove years not of interest:
y.span = c(2000, 2001, 2004, 2006, 2008, 2011, 2013, 2014)
del.span = c(2000:2015)[which(!2000:2015 %in% y.span)]
for(i in del.span){
  del.D <- list.files(path="./tiledD/", pattern=glob2rx(paste0("MOD11A2.A", i,"???.*.005.*.tif$")), full.names=TRUE)
  unlink(del.D)
  del.N <- list.files(path="./tiledN/", pattern=glob2rx(paste0("MOD11A2.A", i,"???.*.005.*.tif$")), full.names=TRUE)
  unlink(del.N)
}

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

aggr_tile <- function(lst, outm, outs){
  if(!file.exists(outm)){
    ## average per block:
    s <- raster::stack(lst)
    smean <- raster::calc(s, fun=mean, na.rm=TRUE)
    ssd <- raster::calc(s, fun=sd, na.rm=TRUE)*10
    r1 <- raster::writeRaster(smean, filename=outm, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
    r2 <- raster::writeRaster(ssd, filename=outs, datatype="INT2U", options=c("COMPRESS=DEFLATE"))
  }
}

## Generate average per month (per year):
aggr_parallel <- function(i){
  lstD <- list.files(path="./tiledD/", pattern=glob2rx(paste0("MOD11A2.*.", modis_grid$TILEh[i], ".005.*.tif$")), full.names=TRUE)
  lstN <- list.files(path="./tiledN/", pattern=glob2rx(paste0("MOD11A2.*.", modis_grid$TILEh[i], ".005.*.tif$")), full.names=TRUE)
  if(length(lstD)>0){
    for(k in 1:length(m.lst)){  
      try( lstD.s <- lstD[unlist(sapply(dirs[which(m.lst[k] == dsel)], function(x){grep(basename(lstD), pattern=glob2rx(paste0("*.A2???",x,".*.*.*.tif$")))}))] )
      if(!class(.Last.value)[1]=="try-error"){
        try( aggr_tile(lstD.s, outm=paste0("./Mtiled/LSTD_M_", m.lst[k], "_", modis_grid$TILEh[i], ".tif"), outs=paste0("./Mtiled/LSTD_SD_", m.lst[k], "_", modis_grid$TILEh[i], ".tif")) )
      }
    }
  }
  if(length(lstN)>0){
    for(k in 1:length(m.lst)){
      try( lstN.s <- lstN[unlist(sapply(dirs[which(m.lst[k] == dsel)], function(x){grep(basename(lstN), pattern=glob2rx(paste0("*.A2???",x,".*.*.*.tif$")))}))] )
      if(!class(.Last.value)[1]=="try-error"){
        try( aggr_tile(lstN.s, outm=paste0("./Mtiled/LSTN_M_", m.lst[k], "_", modis_grid$TILEh[i], ".tif"), outs=paste0("./Mtiled/LSTN_SD_", m.lst[k], "_", modis_grid$TILEh[i], ".tif")) )
      }
    }
  }
}

sfInit(parallel=TRUE, cpus=40)
sfExport("modis_grid", "aggr_tile", "m.lst", "dsel", "dirs")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(sp)
#t <- sfLapply(which(modis_grid$TILEh %in% c("h12v13","h13v13")), aggr_parallel)
t <- sfLapply(1:length(modis_grid$TILEh), aggr_parallel)
sfStop()
## 15,216 tiles

## Mosaic:
M.lstD <- list.files(path="./Mtiled/", pattern="LSTD_M", full.names = TRUE)
M.lstN <- list.files(path="./Mtiled/", pattern="LSTN_M", full.names = TRUE)
SD.lstD <- list.files(path="./Mtiled/", pattern="LSTD_SD", full.names = TRUE)
SD.lstN <- list.files(path="./Mtiled/", pattern="LSTN_SD", full.names = TRUE)
## 3804 / 12 tiles

make_mosaic <- function(j, mon, per, typ){
  outm <- paste0("LST", per[j], "_", typ[j], "_", mon[j], "_1km.tif")
  if(!file.exists(paste0(outm,'.gz'))){
    lst <- get(paste0(typ[j],".lst",per[j]))
    lst.s <- lst[grep(pattern=mon[j], lst)]
    vrt <- paste0(typ[j], '_liste', mon[j], '_', per[j], '.vrt')
    txt <- paste0(typ[j], '_liste', mon[j], '_', per[j], '.txt')
    cat(lst.s, sep="\n", file=txt)
    system(paste0(gdalbuildvrt, ' -input_file_list ', txt,' ', vrt))
    system(paste0(gdalwarp, ' ', vrt, ' ', outm, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"bilinear\" -tr 0.008333333 0.008333333 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))
    #gzip(outm)
    #unlink(vrt)
    #unlink(txt)
  }
}

mon <- rep(m.lst, 4)
per <- as.vector(sapply(c("D","D","N","N"), rep, 12))
typ <- as.vector(sapply(c("M","SD","M","SD"), rep, 12))

sfInit(parallel=TRUE, cpus=38)
sfExport("make_mosaic", "M.lstD", "M.lstN", "SD.lstD", "SD.lstN", "m.lst", "gdalbuildvrt", "gdalwarp", "mon", "per", "typ")
sfLibrary(R.utils)
sfLibrary(sp)
t <- sfLapply(1:48, make_mosaic, mon=mon, per=per, typ=typ)
sfStop()

## Mean annual temperature:
Ml <- paste0("LSTD_M_", m.lst, "_1km.tif")
MNl <- paste0("LSTN_M_", m.lst, "_1km.tif")
#sapply(paste0(Ml,".gz"), gunzip)
#sapply(paste0(MNl,".gz"), gunzip)
sD <- raster::stack(Ml)
sN <- raster::stack(MNl)
## define functions:
meanf <- function(x){calc(x, mean, na.rm=TRUE)}
sdf <- function(x){calc(x, sd, na.rm=TRUE)*10}
## run in parallel:
beginCluster()
r1 <- clusterR(sD, fun=meanf, filename="TMDMOD3.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
## takes ca 15 mins
r2 <- clusterR(sN, fun=meanf, filename="TMNMOD3.tif", datatype="INT2S", options=c("COMPRESS=DEFLATE"))
endCluster()

## Mean 2-month temp:
m2.lst <- c("JanFeb","MarApr","MayJun","JulAug","SepOct","NovDec")
for(i in 1:length(m2.lst)){
  sD <- raster::stack(c(paste0("LSTD_M_", substr(m2.lst[i],1,3), "_1km.tif"), paste0("LSTD_M_", substr(m2.lst[i],4,6), "_1km.tif")))
  fname <- paste0("LSTD_M_", m2.lst[i], "_1km.tif")
  if(!file.exists(fname)){
    beginCluster()
    r1 <- clusterR(sD, fun=meanf, filename=fname, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
    endCluster()
  }
}


## TH: Night time SD images have many missing mixels
## Solution: filter out missing values in all images using land mask
system(paste(gdal_translate, 'TMNMOD3.tif oceanMask1km.sdat -of \"SAGA\"'))
system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"oceanMask1km.sgrd\" -FORMULA=\"ifelse(g1>0,0,1)\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=1 -RESULT=\"oceanMask1km_B.sgrd\"'))
system(paste0('/usr/local/bin/saga_cmd -c=40 grid_calculus 1 -GRIDS=\"oceanMask1km.sgrd\" -FORMULA=\"ifelse(g1>0,1,0)\" -INTERPOLATION=0 -USE_NODATA=1 -TYPE=1 -RESULT=\"landMask1km_B.sgrd\"'))
GDALinfo("oceanMask1km_B.sdat")
unlink("oceanMask1km_B.7z")
system(paste0("7za a -t7z oceanMask1km_B.7z oceanMask1km_B.* -mmt=40"))
unlink("landMask1km_B.7z")
system(paste0("7za a -t7z landMask1km_B.7z landMask1km_B.* -mmt=40"))

SNl <- paste0("LSTN_SD_", m.lst, "_1km.tif")
na.count.mask <- sum(getValues(is.na(raster("landMask1km_B.sdat")))) ## TAKES 3-4 mins!

for(i in 1:length(SNl)){
  out <- gsub("1km.tif", "1kmf.tif", SNl[i])
  if(!file.exists(out)){
    na.count <- sum(getValues(is.na(raster(SNl[i]))))
    if(na.count>na.count.mask){
      tmp <- set.file.extension(tempfile(), ".sdat")
      tmp2 <- set.file.extension(tempfile(), ".sgrd")
      system(paste0(gdal_translate, ' ', SNl[i], ' ', tmp, ' -of \"SAGA\"'))
      system(paste0('/usr/local/bin/saga_cmd -c=40 grid_tools 29 -INPUT ', set.file.extension(tmp, ".sgrd"), ' -MASK=\"landMask1km_B.sgrd\" -INTERPOLATION=1 -GROW=2 -PYRAMIDS=1 -START=0 -RESULT=\"', tmp2, '\"'))
      system(paste0(gdal_translate, ' ', set.file.extension(tmp2, ".sdat"), ' ', out, ' -ot \"Int16\" -a_srs \"+proj=longlat +datum=WGS84\" -co \"COMPRESS=DEFLATE\"')) ## -a_nodata 65535
      unlink(tmp)
      unlink(set.file.extension(tmp2, ".sdat"))
    }
  }
}


## Distance to the sea (global)
## reproject to Equi7 system (1 km):
load("/data/EcoTapestry/equi7t3.rda")
sfInit(parallel=TRUE, cpus=7)
sfLibrary(sp)
sfExport("equi7t3", "gdalwarp")
x <- sfLapply(1:length(equi7t3), function(x){system(paste0(gdalwarp, ' oceanMask1km_B.sdat OCM_', names(equi7t3)[x], '_2km.sdat -t_srs \"', proj4string(equi7t3[[x]]), '\" -tr 2000 2000 -r \"near\" -of \"SAGA\" -dstnodata 255 -te ', paste(as.vector(bbox(equi7t3[[x]])), collapse=" "), ' -ot \"Byte\"')) })

## Buffer distances (Euclidean / in meters):
sfInit(parallel=TRUE, cpus=7)
sfExport("equi7t3", "gdalwarp")
x <- sfLapply(1:length(equi7t3), function(x){ if(!file.exists(paste0('DSTOCE3_', names(equi7t3)[x], '.sgrd,'))){ system(paste0('/usr/local/bin/saga_cmd -c=1 grid_tools 26 -FEATURES=\"OCM_', names(equi7t3)[x], '_2km.sgrd\" -DISTANCE=\"DSTOCM_', names(equi7t3)[x], '.sgrd\"')) }})
sfStop()

