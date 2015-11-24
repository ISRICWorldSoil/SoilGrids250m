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
## 46 dirs
days <- as.numeric(format(seq(ISOdate(2015,1,1), ISOdate(2015,12,31), by="month"), "%j"))-1
m.lst <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dsel <- cut(as.numeric(dirs), breaks=c(days,366), labels = m.lst)
modis.prj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
xs1 <- list.files(path="./tilednir", pattern=glob2rx("*.A2014257.*.*.*.tif$"))
TILEh1 <- levels(as.factor(sapply(xs, function(x){strsplit(x, "\\.")[[1]][3]})))
## 326 x 12 MODIS tiles

aggr_tile <- function(lst, outm){
  if(!file.exists(outm)){
    s <- raster::stack(lst)
    smean <- raster::calc(s, fun=mean, na.rm=TRUE)*10000
    r1 <- raster::writeRaster(smean, filename=outm, datatype="INT2S", options=c("COMPRESS=DEFLATE", overwrite=TRUE))
    gc()
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
aggr_parallel(i=421, channel="mir", band=7)
aggr_parallel(i=597, channel="nir", band=4)
#aggr_parallel(i=195, channel="nir", band=4)

## TAKES CA 2 hrs
sfInit(parallel=TRUE, cpus=40)
sfExport("modis_grid", "aggr_tile", "m.lst", "dsel", "dirs")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(sp)
t <- sfLapply(1:length(modis_grid$TILEh), aggr_parallel, channel="nir", band=4)
sfStop()

sfInit(parallel=TRUE, cpus=40)
sfExport("modis_grid", "aggr_tile", "m.lst", "dsel", "dirs")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(sp)
t <- sfLapply(1:length(modis_grid$TILEh), aggr_parallel, channel="mir", band=7)
sfStop()

## list of output tiles:
M4.lstD <- list.files(path="./Mtilednir", pattern=glob2rx("NRB4_M_*_*.tif$"), full.names = TRUE)
TILEh2 <- levels(as.factor(sapply(M4.lstD, function(x){strsplit(basename(x), "_")[[1]][4]})))
#v01.rm <- M4.lstD[grep(pattern="v01", M4.lstD)]
M7.lstD <- list.files(path="./Mtiledmir", pattern=glob2rx("NRB7_M_*_*.tif$"), full.names = TRUE)
## 3864 tiles x 2


## TH: Some months miss many pixels, especially Jan and Feb for latitudes above 63 degree and Jun Jul and Aug for tropics.
## To prepare the land mask at 500 m we use:
t.lst <- M4.lstD[grep(pattern="_M_Jan", M4.lstD)]
## 322
mod.land <- lapply(t.lst, function(x){c(as.vector(extent(raster(x))), dim(raster(x)))})
modis_grid_land500m <- as.data.frame(do.call(rbind, mod.land))
names(modis_grid_land500m) <- c("xmin","xmax","ymin","ymax","nrow","ncol")
modis_grid_land500m$TILEh <- sapply(t.lst, function(x){strsplit(strsplit(x, "_")[[1]][4], "\\.")[[1]][1]})
summary(!modis_grid_land500m$nrow==2400)
save(modis_grid_land500m, file="modis_grid_land500m.rda")

## Replace missing values using the neighbouring months
## (this function can be run multiple times over the same tiles)
filter_NA <- function(i, mask.path="/data/modis_grid"){
  c.lst <- rep(m.lst, 3)
  m <- raster(paste0(mask.path, "/LMK_", i, "_500m.tif"))
  ## run only for tiles with land mask
  if(!sum(getValues(is.na(m)))==length(m)){
    tgrid <- raster(paste0("./Mtilednir/NRB4_M_", m.lst[1],"_", i, ".tif"))
    ## some tiles got 1 pixel extra and need to be re-aligned
    if(!length(tgrid)==length(m)){
      ## if necessary resample:
      m <- raster::resample(m, tgrid, method='ngb')      
    }
    m <- as(as(m, "SpatialGridDataFrame"), "SpatialPixelsDataFrame") 
    for(j in 1:length(m.lst)){  
      infile1 <- paste0("./Mtilednir/NRB4_M_", m.lst[j],"_", i, ".tif")
      infile2 <- paste0("./Mtiledmir/NRB7_M_", m.lst[j],"_", i, ".tif")
      nlst <- c(c.lst[j+12-2], c.lst[j+12-1], c.lst[j+12+1], c.lst[j+12+2])
      ## read target layer:
      m$d1 <- readGDAL(infile1, silent=TRUE)$band1[m@grid.index]
      ## check if there are missing pixels, if yes replace with temporal neighbours:
      if(any(is.na(m$d1))){
        try( nb1 <- readGDAL(paste0("./Mtilednir/NRB4_M_", nlst[1], "_", i, ".tif"))$band1[m@grid.index] )
        if(!exists(nb1)){ nb1b = rep(NA, length(m$d1)) }
        try( nb2 <- readGDAL(paste0("./Mtilednir/NRB4_M_", nlst[2], "_", i, ".tif"))$band1[m@grid.index] )
        if(!exists(nb2)){ nb2b = rep(NA, length(m$d1)) }
        try( nb3 <- readGDAL(paste0("./Mtilednir/NRB4_M_", nlst[3], "_", i, ".tif"))$band1[m@grid.index] )
        if(!exists(nb3)){ nb3b = rep(NA, length(m$d1)) }
        try( nb4 <- readGDAL(paste0("./Mtilednir/NRB4_M_", nlst[4], "_", i, ".tif"))$band1[m@grid.index] )
        if(!exists(nb4)){ nb4b = rep(NA, length(m$d1)) }
        m$df1 <- rowMeans(cbind(nb1, nb2, nb2, nb3, nb3, nb4), na.rm=TRUE)
        m$df1 <- ifelse(is.na(m$d1), m$df1, m$d1)
        writeGDAL(m["df1"], infile1, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      }
      m$d2 <- readGDAL(infile2, silent=TRUE)$band1[m@grid.index]
      if(any(is.na(m$d2))){
        try( nb1b <- readGDAL(paste0("./Mtiledmir/NRB7_M_", nlst[1], "_", i, ".tif"))$band1[m@grid.index] )
        if(!exists(nb1b)){ nb1b = rep(NA, length(m$d2)) }
        try( nb2b <- readGDAL(paste0("./Mtiledmir/NRB7_M_", nlst[2], "_", i, ".tif"))$band1[m@grid.index] )
        if(!exists(nb2b)){ nb2b = rep(NA, length(m$d2)) }
        try( nb3b <- readGDAL(paste0("./Mtiledmir/NRB7_M_", nlst[3], "_", i, ".tif"))$band1[m@grid.index] )
        if(!exists(nb3b)){ nb3b = rep(NA, length(m$d2)) }
        try( nb4b <- readGDAL(paste0("./Mtiledmir/NRB7_M_", nlst[4], "_", i, ".tif"))$band1[m@grid.index] )
        if(!exists(nb4b)){ nb4b = rep(NA, length(m$d2)) }
        m$df2 <- rowMeans(cbind(nb1b, nb2b, nb2b, nb3b, nb3b, nb4b), na.rm=TRUE)
        m$df2 <- ifelse(is.na(m$d2), m$df2, m$d2)
        writeGDAL(m["df2"], infile2, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      }
      gc()
    }
  }
}

## test it:
filter_NA(i="h20v01")

sfInit(parallel=TRUE, cpus=40)
sfExport("filter_NA", "m.lst", "modis_grid_land500m")
sfLibrary(rgdal)
sfLibrary(raster)
t <- sfClusterApplyLB(modis_grid$TILEh, function(i){try( filter_NA(i) )})
sfStop()

## Mosaic:
make_mosaic <- function(j, band){
  outm <- paste0("NRB", band, "_M_", m.lst[j], "_500m.tif")
  if(!file.exists(paste0(outm,'.gz'))){
    lst <- get(paste0("M", band, ".lstD"))
    lst.s <- lst[grep(pattern=m.lst[j], lst)]
    vrt <- paste0('M_liste', m.lst[j], '_', band, '.vrt')
    txt <- paste0('M_liste', m.lst[j], '_', band, '.txt')
    cat(lst.s, sep="\n", file=txt)
    system(paste0(gdalbuildvrt, ' -input_file_list ', txt,' ', vrt))
    #system(paste0(gdalwarp, ' ', vrt, ' ', outm, ' -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Int16\" -dstnodata \"-32768\" -tr 0.004166667 0.004166667 -te -180 -90 180 90 -multi'))
    #gzip(outm)
    #unlink(vrt)
    #unlink(txt)
  }
}

## test:
make_mosaic(j=2, band=7)

bands <- c(rep(4,12), rep(7,12))
mon <- rep(1:length(m.lst), 2)
for(x in 1:length(mon)){
  make_mosaic(mon[x], band=bands[x])
}

# sfInit(parallel=TRUE, cpus=6)
# sfExport("make_mosaic", "M4.lstD", "M7.lstD", "m.lst", "gdalbuildvrt", "gdalwarp", "mon", "bands")
# sfLibrary(R.utils)
# sfLibrary(sp)
# t <- sfLapply(1:length(mon), function(x){ make_mosaic(mon[x], band=bands[x])})
# sfStop()


## Create world map of shifting sand / bedrock:
deserts.sin <- function(j, grid=modis_grid_land500m, t_srs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m+no_defs"){
  nfile <- paste0("extra/desertPR_", grid$TILEh[j], ".tif")
  nfile2 <- paste0("extra/barerockPR_", grid$TILEh[j], ".tif")
  te <- as.vector(grid[j,c(1,3,2,4)])
  if(!file.exists(nfile)|!file.exists(nfile2)){  
    dfile <- paste0("extra/LSTD_", grid$TILEh[j], ".tif")
    if(!file.exists(dfile)){
      system(paste0(gdalwarp, ' /data/MOD10A2/LSTD_M_annual_500m.sdat ', dfile, ' -t_srs \"', t_srs, '\" -ts 2400 2400 -te ', paste(te, collapse=" "), ' -co \"COMPRESS=DEFLATE\"'))  ## -r \"near\" -ot \"Byte\"
    }
    kfile <- paste0("extra/SLP_", grid$TILEh[j], ".tif")
    if(!file.exists(kfile)){
      system(paste0(gdalwarp, ' /data/MDEM/SLPSRM3a.tif ', kfile, ' -t_srs \"', t_srs, '\" -ts 2400 2400 -te ', paste(te, collapse=" "), ' -co \"COMPRESS=DEFLATE\"'))
    }
    ## probability of shifting sand:
    lfiles <- M7.lstD[grep(pattern=grid$TILEh[j], M7.lstD)]
    try(  m <- as(raster::stack(c(dfile, kfile, lfiles)), "SpatialGridDataFrame") )
    if(!class(.Last.value)[1]=="try-error"){
      MIR <- rowMeans(m@data[,-c(1:2)], na.rm=TRUE)
      sel <- MIR>4050 & m@data[,1] > 298 & m@data[,2] < 0.15
      ## hist(m@data[m@data[,1]>0,1])
      ## this needs to be finetuned (or use some ground data?)
      sel2 <- MIR>1300 & m@data[,2] > 0.32
      if(sum(sel, na.rm=TRUE)>0){
        if(!file.exists(nfile)){
          m$d <- ifelse(sel, 1, NA) 
          writeGDAL(m["d"], nfile, type="Byte", mvFlag=0, options="COMPRESS=DEFLATE")
        }
      }
      if(sum(sel2, na.rm=TRUE)>0){
        if(!file.exists(nfile2)){
          m$r <- ifelse(sel2, 1, NA) 
          writeGDAL(m["r"], nfile2, type="Byte", mvFlag=0, options="COMPRESS=DEFLATE")
        }
      }
    }
  }
}
sfInit(parallel=TRUE, cpus=40)
sfLibrary(rgdal)
sfLibrary(raster)
sfExport("modis_grid_land500m", "deserts.sin", "gdalwarp", "M7.lstD")
x <- sfLapply(1:nrow(modis_grid_land500m), function(j){try(deserts.sin(j))})
sfStop()

pr2.lst <- list.files(path="extra/", pattern=glob2rx("barerockPR_*.tif$"), full.names=TRUE)
unlink("my_liste_PR2.txt")
cat(pr2.lst, sep="\n", file="my_liste_PR2.txt")
system(paste(gdalbuildvrt, '-input_file_list my_liste_PR2.txt barerockPR.vrt'))
unlink("barerockPR_sin.tif")
system(paste(gdalwarp, 'barerockPR.vrt barerockPR_sin.tif -tr 500 500 -ot \"Byte\" -dstnodata 0 -co \"COMPRESS=DEFLATE\"'))

pr.lst <- list.files(path="extra/", pattern=glob2rx("desertPR_*.tif$"), full.names=TRUE)
unlink("my_liste_PR.txt")
cat(pr.lst, sep="\n", file="my_liste_PR.txt")
system(paste(gdalbuildvrt, '-input_file_list my_liste_PR.txt desertPR.vrt'))
unlink("desertPR_sin.tif")
system(paste(gdalwarp, 'desertPR.vrt desertPR_sin.tif -tr 500 500 -ot \"Byte\" -dstnodata 0 -co \"COMPRESS=DEFLATE\"'))

