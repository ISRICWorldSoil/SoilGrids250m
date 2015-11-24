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
## only 296 cover land
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

## filter tiles with missing values for landmask:
filter_NA <- function(i, mask.path="/data/modis_grid"){
  ## if necessary resample:
  tgrid <- raster("./Mtiled/EVI_M_", m.lst[1],"_", i, ".tif")
  m <- raster(paste0(mask.path, "/LMK_", i, ".tif"))
  if(!sum(getValues(is.na(m)))==length(m)){
    ## some tiles got 1 pixel extra and need to be re-aligned
    if(!length(tgrid)==length(m)){
      m <- raster::resample(m, tgrid, method='ngb')      
    }
    m <- as(as(m, "SpatialGridDataFrame"), "SpatialPixelsDataFrame") ## takes 1-2 mins
    for(j in 1:length(m.lst)){
      infile1 <- paste0("./Mtiled/EVI_M_", m.lst[j],"_", i, ".tif")
      infile2 <- paste0("./Mtiled/EVI_SD_", m.lst[j],"_", i, ".tif")
      n1 <- ifelse(j==1, m.lst[length(m.lst)], m.lst[j-1])
      n2 <- ifelse(j==length(m.lst), m.lst[1], m.lst[j+1])
      m$d1 <- readGDAL(infile1, silent=TRUE)$band1[m@grid.index]
      ## check if there are missing pixels:
      if(any(is.na(m$d1))){
        m$df1 <- rowMeans(cbind(readGDAL(paste0("./Mtiled/EVI_M_", n1, "_", i, ".tif"))$band1[m@grid.index], readGDAL(paste0("./Mtiled/EVI_M_", n2, "_", i, ".tif"))$band1[m@grid.index], silent=TRUE), na.rm=TRUE)
        m$df1 <- ifelse(is.na(m$d1), m$df1, m$d1)
        writeGDAL(m["df1"], infile1, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      }
      m$d2 <- readGDAL(infile2, silent=TRUE)$band1[m@grid.index]
      if(any(is.na(m$d2))){
        m$df2 <- rowMeans(cbind(readGDAL(paste0("./Mtiled/EVI_SD_", n1, "_", i, ".tif"))$band1[m@grid.index], readGDAL(paste0("./Mtiled/EVI_SD_", n2, "_", i, ".tif"))$band1[m@grid.index], silent=TRUE), na.rm=TRUE)
        m$df2 <- ifelse(is.na(m$d2), m$df2, m$d2)
        writeGDAL(m["df2"], infile2, type="Int16", mvFlag=-32768, options="COMPRESS=DEFLATE")
      }
      gc()    
    }
  }
}

sfInit(parallel=TRUE, cpus=40)
sfExport("filter_NA", "m.lst", "modis_grid_land")
sfLibrary(rgdal)
sfLibrary(raster)
t <- sfLapply(modis_grid_land$TILEh, function(i){try( filter_NA(i) )})
sfStop()

## Mosaic:
M.lstD <- list.files(path="./Mtiled", pattern=glob2rx("EVI_M_*_*.tif$"), full.names = TRUE)
SD.lstD <- list.files(path="./Mtiled", pattern=glob2rx("EVI_SD_*_*.tif$"), full.names = TRUE)
## 3552 tiles
JanFeb.lst <- M.lstD[grep(pattern="_M_JanFeb", M.lstD)]
## get bounding box and resolution
mod.land <- lapply(JanFeb.lst, function(x){c(as.vector(extent(raster(x))), res(raster(x)))})
#mod.land <- unique(sapply(basename(M.lstD[grep(pattern="_M_JanFeb", M.lstD)]), function(x){strsplit(strsplit(x, "_")[[1]][4], "\\.")[[1]][1]}))
#modis_grid_land <- modis_grid[which(modis_grid$TILEh %in% mod.land),]
modis_grid_land <- as.data.frame(do.call(rbind, mod.land))
names(modis_grid_land) <- c("xmin","xmax","ymin","ymax","xres","yres")
modis_grid_land$TILEh <- sapply(JanFeb.lst, function(x){strsplit(strsplit(x, "_")[[1]][4], "\\.")[[1]][1]})
save(modis_grid_land, file="modis_grid_land.rda")


make_mosaic <- function(j, mon, typ, res=c("250m","1km")){
  outm <- paste0("EVI", "_", typ[j], "_", mon[j], "_", res[1], ".tif")
  lst <- get(paste0(typ[j],".lstD"))
  lst.s <- lst[grep(pattern=mon[j], lst)]
  vrt <- paste0(typ[j], '_liste', mon[j], '.vrt')
  txt <- paste0(typ[j], '_liste', mon[j], '.txt')
  unlink(txt)
  cat(lst.s, sep="\n", file=txt)
  system(paste0(gdalbuildvrt, ' -input_file_list ', txt,' ', vrt))
  if(any(res=="250m")){
    system(paste0(gdalwarp, ' ', vrt, ' ', outm, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"bilinear\" -ot \"Int16\" -dstnodata \"-32768\" -tr 0.002083333 0.002083333 -te -180 -90 180 90 -co BIGTIFF=YES -wm 4000')) ## -co \"COMPRESS=DEFLATE\" 
    gzip(outm)
  }
  if(any(res=="1km")){
    system(paste0(gdalwarp, ' ', vrt, ' ', gsub("250m", "1km", outm), ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"bilinear\" -ot \"Int16\" -dstnodata \"-32768\" -tr 0.008333333 0.008333333 -te -180 -90 180 90 -co \"COMPRESS=DEFLATE\"'))
  }
  #unlink(vrt)
  #unlink(txt)
}

mon <- rep(m.lst, 2)
typ <- as.vector(sapply(c("M","SD"), rep, length(m.lst)))
## in series:
for(j in 1:length(mon)){
  make_mosaic(j, mon, typ, res="250m")
}

sfInit(parallel=TRUE, cpus=6) ## too many processes lead to error ('No space left on device')
sfExport("make_mosaic", "M.lstD", "SD.lstD", "m.lst", "gdalbuildvrt", "gdalwarp", "mon", "typ")
sfLibrary(R.utils)
sfLibrary(rgdal)
t <- sfLapply(1:length(mon), make_mosaic, mon=mon, typ=typ, res="1km")
sfStop()
