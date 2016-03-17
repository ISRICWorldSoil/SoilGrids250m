# title         : soilmask_250m.R
# purpose       : Preparation of the global soil mask;
# reference     : SoilGrids250m
# producer      : Prepared by T. Hengl (tom.hengl@wur.nl)
# address       : In Wageningen, NL.
# inputs        : GLC30 zipped files download from http://www.globallandcover.com/GLC30Download/download.aspx
# outputs       : 250 m resolution imagery;
# remarks 1     : 853 blocks;

library(RSAGA)
library(gdalUtils)
library(rgdal)
library(R.utils)
library(utils)
library(snowfall)
library(raster)
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate = "gdal_translate"
  gdalwarp = "gdalwarp"
}
load("E:\\EQUI7\\equi7t3.rda")
load("E:\\EQUI7\\equi7t1.rda")

classes <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
classes.t <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10")
zip.lst <- list.files(path="X:\\GlobCover30\\2010_original", pattern=glob2rx("*.zip$"), full.names=TRUE)
str(zip.lst)


resample.glc30 <- function(i, classes=classes){
  try( lst <- unzip(zip.lst[i], list = TRUE) )
  if(!class(.Last.value)[1]=="try-error"){
    inname <- lst$Name[grep(lst$Name, pattern=glob2rx("*.tif$"))]   
    name30m <- gsub(".tif", paste0("_30m.tif"), inname)
    if(!file.exists(name30m)){
      if(nchar(inname)>0){
        unzip(zip.lst[i], inname)
        ## resample to local coordinates:
        system(paste(gdalwarp, inname, gsub(".tif", ".sdat", inname), '-of \"SAGA\" -r \"near\" -srcnodata 0 -ot \"Byte\" -dstnodata 255 -overwrite'))
        system(paste(gdal_translate, gsub(".tif", ".sdat", inname), name30m, '-co \"COMPRESS=DEFLATE\" -ot \"Byte\"'))
        ## mask values per each class:
        for(k in 1:length(classes)){
          outn <- gsub(".tif", paste0("_", classes[k] ,"_250m.tif"), inname)
          rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS=gsub(".tif", ".sgrd", inname), RESULT=gsub(".tif", paste0("_", classes[k] ,".sgrd"), inname), FORMULA=paste("ifelse(a=",classes[k],",100,0)", sep="")), show.output.on.console = FALSE)
          ## resample:
          system(paste(gdalwarp, gsub(".tif", paste0("_", classes[k] ,".sdat"), inname), outn, '-r \"average\" -ot \"Byte\" -tr 250 250 -co \"COMPRESS=DEFLATE\" -dstnodata 255 -overwrite'))
          ## clean up:
          unlink(gsub(".tif", paste0("_", classes[k], ".sgrd"), inname))
          unlink(gsub(".tif", paste0("_", classes[k], ".mgrd"), inname))
          unlink(gsub(".tif", paste0("_", classes[k], ".sdat"), inname))
          unlink(gsub(".tif", paste0("_", classes[k], ".prj"), inname))
        }
        unlink(gsub(".tif", ".sdat", inname))
        unlink(gsub(".tif", ".sgrd", inname))
        unlink(gsub(".tif", ".prj", inname))
        unlink(gsub(".tif", ".sdat.aux.xml", inname))
        unlink(inname)
      }
    }
  }
}

lapply(1:length(zip.lst), resample.glc30, classes=classes)

#sfInit(parallel=TRUE, cpus=6)
#sfLibrary(RSAGA)
#sfLibrary(utils)
#sfLibrary(gdalUtils)
#sfExport("zip.lst", "classes", "csize", "af.csy")
#t <- sfLapply(1:length(zip.lst), resample.glc30, classes=classes)
#sfStop()


## Global mosaic at 1/480 resolution (about 250 m) is 172,800 cols x 86,400 rows

## Create mosaick:
outdir <- "G:\\SoilGrids250m\\zipped"

## list all geotifs, create mosaicks and compress to output dir:
tifout.lst <- paste0("*_", (1:10)*10, "_250m.tif$")
for(j in 1:length(tifout.lst)){
  tmp.lst <- list.files(pattern=glob2rx(tifout.lst[j]), recursive = TRUE)
  ## 852 tiles per layer!
  ## reproject to the same coord system:
  sfInit(parallel=TRUE, cpus=6)
  sfExport("tmp.lst", "gdalwarp")
  t <- sfLapply(tmp.lst, function(x){ if(!file.exists(paste0("tmp/", gsub("_250m.tif", "_ll.tif", basename(x))))){ system(paste0(gdalwarp, ' ', x, ' -t_srs \"+proj=longlat +datum=WGS84\" ', paste0("tmp/", gsub("_250m.tif", "_ll.tif", basename(x))), ' -tr 0.002083333 0.002083333 -ot \"Byte\" -dstnodata 255')) } })
  sfStop() ## end of parallel processing...
}

## Derive soil mask using a simple formula:
GLC_90 <- list.files(path="./tmp", pattern=glob2rx("*_*_*_90_ll.tif$"), full.names=TRUE) 
GLC_60 <- list.files(path="./tmp", pattern=glob2rx("*_*_*_60_ll.tif$"), full.names=TRUE)
GLC_100 <- list.files(path="./tmp", pattern=glob2rx("*_*_*_100_ll.tif$"), full.names=TRUE)
GLC.lst <- data.frame(GLC_90, GLC_60, GLC_100)
str(GLC.lst)
## 3 x 852 tiles

makeSmask <- function(GLC_90, GLC_60, GLC_100){
  outName <- paste0("tmp/", gsub("_90_", "_SMK_", basename(paste(GLC_90))))
  if(!file.exists(outName)){
    gc()
    try( m <- stack(c(paste(GLC_90), paste(GLC_60), paste(GLC_100))) )
    if(!class(.Last.value)[1]=="try-error"){
      m <- as(m, "SpatialGridDataFrame")
      names(m) <- c("GLC_90", "GLC_60", "GLC_100")
      try( m$mask <- ifelse(m$GLC_60>90|m$GLC_90>90|m$GLC_100>90, 0, 1) )
      try( writeGDAL(m["mask"], outName, type="Byte", mvFlag=0) ) # , options="COMPRESS=DEFLATE"
      gc()
    }
  }
}
## test:
#makeSmask(GLC_90=GLC.lst[138,1], GLC_60=GLC.lst[138,2], GLC_100=GLC.lst[138,3])

sfInit(parallel=TRUE, cpus=6)
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(raster)
sfExport("GLC.lst", "makeSmask")
x <- sfLapply(1:nrow(GLC.lst), function(x){makeSmask(GLC_90=GLC.lst[x,1], GLC_60=GLC.lst[x,2], GLC_100=GLC.lst[x,3])} )
sfStop()

## Land mask - includes deserts!
makeLmask <- function(GLC_60, GLC_100){
  outName <- paste0("tmp/", gsub("_60_", "_LMK_", basename(paste(GLC_60))))
  if(!file.exists(outName)){
    gc()
    try( m <- stack(c(paste(GLC_60), paste(GLC_100))) )
    if(!class(.Last.value)[1]=="try-error"){
      m <- as(m, "SpatialGridDataFrame")
      names(m) <- c("GLC_60", "GLC_100")
      try( m$mask <- ifelse(m$GLC_60>90|m$GLC_100>90, 0, 1) )
      try( writeGDAL(m["mask"], outName, type="Byte", mvFlag=0) ) # , options="COMPRESS=DEFLATE"
      gc()
    }
  }
}

sfInit(parallel=TRUE, cpus=6)
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(raster)
sfExport("GLC.lst", "makeLmask")
x <- sfLapply(1:nrow(GLC.lst), function(x){makeLmask(GLC_60=GLC.lst[x,2], GLC_100=GLC.lst[x,3])} )
sfStop()

## per continent:
smk.lst <- list.files(path="tmp/", pattern=glob2rx("*_SMK_ll.tif$"), full.names=TRUE)
unlink("my_liste.txt")
cat(smk.lst, sep="\n", file="my_liste.txt")
gdalbuildvrt(input_file_list="my_liste.txt", output.vrt="smask.vrt")

tile.tif <- function(t, t_srs = proj4string(t)){
  for(j in 1:nrow(t)){
    nfile <- paste0("stiled/SMK_", strsplit(paste(t@data[j,"SHORTNAME"]), " ")[[1]][2], "_", t@data[j,"TILE"], ".tif")
    te <- as.vector(bbox(t[j,]))
    system(paste0(gdalwarp, ' smask.vrt ', nfile, ' -t_srs \"', t_srs, '\" -tr 250 250 -r \"near\" -srcnodata 0 -dstnodata 0 -te ', paste(te, collapse=" "), ' -ot \"Byte\"'))
  }
}
sfInit(parallel=TRUE, cpus=6)
sfLibrary(sp)
sfExport("equi7t3", "tile.tif", "gdalwarp")
x <- sfLapply(equi7t3, tile.tif)
sfStop() ## end of parallel processing

## Land mask in MODIS sinusoidal projection
load("../MOD13Q1/modis_grid_land.rda")
str(modis_grid_land)
tile.sin <- function(j, grid=modis_grid_land, t_srs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m+no_defs"){
  nfile <- paste0("sinusoidal/LMK_", grid$TILEh[j], ".tif")
  te <- as.vector(grid[j,c(1,3,2,4)])
  if(!file.exists(nfile)){
    system(paste0(gdalwarp, ' lmask.vrt ', nfile, ' -t_srs \"', t_srs, '\" -ts 4800 4800 -r \"near\" -srcnodata 0 -dstnodata 0 -te ', paste(te, collapse=" "), ' -ot \"Byte\" -co \"COMPRESS=DEFLATE\"')) ## 
  }
  dfile <- paste0("sinusoidal/DEM_", grid$TILEh[j], ".tif")
  if(!file.exists(dfile)){
    system(paste0(gdalwarp, ' H:/srtm15_plus/GMTED2010_250m.tif ', dfile, ' -t_srs \"', t_srs, '\" -ts 4800 4800 -te ', paste(te, collapse=" "), ' -co \"COMPRESS=DEFLATE\"'))  ## -r \"near\" -ot \"Byte\"
  }
  kfile <- paste0("sinusoidal/MOD44W_", grid$TILEh[j], ".tif")
  if(!file.exists(kfile)){
    system(paste0(gdalwarp, ' E:/WORLDGRIDS/MOD44W/MOD44W.vrt ', kfile, ' -t_srs \"', t_srs, '\" -ts 4800 4800 -te ', paste(te, collapse=" "), ' -dstnodata 255 -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))
    try( m <- as(stack(c(nfile, dfile, kfile)), "SpatialGridDataFrame") )
    if(!class(.Last.value)[1]=="try-error"){
      ## Filter lines and artifacts in the sea (GLC30):
      sel <- m@data[,1]==1 & ( !m@data[,3]==100 | (!m@data[,2]==0|is.na(m@data[,2])) )
      if(sum(sel, na.rm=TRUE)>0){
        m$d <- ifelse(sel, 1, NA) ## | is.na(m@data[,3])
        writeGDAL(m["d"], nfile, type="Byte", mvFlag=0, options="COMPRESS=DEFLATE")
      }
    }
  }
}
sfInit(parallel=TRUE, cpus=3)
sfLibrary(sp)
sfExport("modis_grid_land", "tile.sin", "gdalwarp")
x <- sfLapply(1:nrow(modis_grid_land), function(j){try(tile.sin(j))})
sfStop() ## end of parallel processing

## Land mask in MODIS sinusoidal projection 500 m resolution
load("../MCD43A4/modis_grid_land500m.rda")
str(modis_grid_land500m)
lst.rem <- list.files(path="./sinusoidal", pattern=glob2rx("LMK_*_500m.tif"), full.names = TRUE)
unlink(lst.rem)

tile.sin500m <- function(j, grid=modis_grid_land500m, t_srs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m+no_defs", replace=TRUE){
  nfile <- paste0("sinusoidal/LMK_", grid$TILEh[j], "_500m.tif")
  te <- as.vector(grid[j,c(1,3,2,4)])
  if(!file.exists(nfile)){
    system(paste0(gdalwarp, ' lmask.vrt ', nfile, ' -t_srs \"', t_srs, '\" -ts 2400 2400 -r \"near\" -srcnodata 0 -dstnodata 0 -te ', paste(te, collapse=" "), ' -ot \"Byte\" -co \"COMPRESS=DEFLATE\"')) ## 
  }
  dfile <- paste0("sinusoidal/DEM_", grid$TILEh[j], "_500m.tif")
  if(!file.exists(dfile)){
    system(paste0(gdalwarp, ' H:/srtm15_plus/GMTED2010_250m.tif ', dfile, ' -t_srs \"', t_srs, '\" -ts 2400 2400 -te ', paste(te, collapse=" "), ' -co \"COMPRESS=DEFLATE\"'))  ## -r \"near\" -ot \"Byte\"
  }
  kfile <- paste0("sinusoidal/MOD44W_", grid$TILEh[j], "_500m.tif")
  if(!file.exists(kfile)){
    system(paste0(gdalwarp, ' E:/WORLDGRIDS/MOD44W/MOD44W.vrt ', kfile, ' -t_srs \"', t_srs, '\" -ts 2400 2400 -te ', paste(te, collapse=" "), ' -dstnodata 255 -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))
  }
  if(replace==TRUE){  
    try( m <- as(stack(c(nfile, dfile, kfile)), "SpatialGridDataFrame") )
    if(!class(.Last.value)[1]=="try-error"){
      ## Filter lines and artifacts in the sea (GLC30):
      sel <- m@data[,1]==1 & ( !m@data[,3]==100 | (!m@data[,2]==0|is.na(m@data[,2])) )
      if(sum(sel, na.rm=TRUE)>0){
        m$d <- ifelse(sel, 1, NA) ## | is.na(m@data[,3])
        writeGDAL(m["d"], nfile, type="Byte", mvFlag=0, options="COMPRESS=DEFLATE")
      }
    }
  }
}
sfInit(parallel=TRUE, cpus=10)
sfLibrary(sp)
sfExport("modis_grid_land500m", "tile.sin500m", "gdalwarp")
x <- sfLapply(1:nrow(modis_grid_land500m), tile.sin500m)
sfStop() ## end of parallel processing

## landmask per continent:
lmk.lst <- list.files(path="tmp/", pattern=glob2rx("*_LMK_ll.tif$"), full.names=TRUE)
unlink("my_liste.txt")
cat(lmk.lst, sep="\n", file="my_liste.txt")
gdalbuildvrt(input_file_list="my_liste.txt", output.vrt="lmask.vrt")

tile2.tif <- function(t, t_srs = proj4string(t)){
  for(j in 1:nrow(t)){
    cn = strsplit(paste(t@data[j,"SHORTNAME"]), " ")[[1]][2]
    te <- as.vector(bbox(t[j,]))
    nfile <- paste0("stiled/LMK_", cn, "_", t@data[j,"TILE"], ".tif")
    if(!file.exists(nfile)){
      system(paste0(gdalwarp, ' lmask.vrt ', nfile, ' -t_srs \"', t_srs, '\" -tr 250 250 -r \"near\" -srcnodata 0 -dstnodata 0 -te ', paste(te, collapse=" "), ' -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))
    }
    dfile <- paste0("stiled/DEM_", cn, "_", t@data[j,"TILE"], ".tif")
    if(!file.exists(dfile)){
      system(paste0(gdalwarp, ' X:/DEMs/MDEM_', cn, '_250m.sdat ', dfile, ' -t_srs \"', t_srs, '\" -tr 250 250 -te ', paste(te, collapse=" "), ' -co \"COMPRESS=DEFLATE\"')) 
    }
    kfile <- paste0("stiled/MOD44W_", cn, "_", t@data[j,"TILE"], ".tif")
    if(!file.exists(kfile)){
      system(paste0(gdalwarp, ' E:/WORLDGRIDS/MOD44W/MOD44W.vrt ', kfile, ' -t_srs \"', t_srs, '\" -tr 250 250 -te ', paste(te, collapse=" "), ' -dstnodata 255 -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))  ## -r \"bilinear\"
      ## Filter lines and artifacts in the sea (GLC30):
      ## needs to be land mask by at least 2 sources:
      if(!cn=="AN"){
        try( m <- as(stack(c(nfile, dfile, kfile)), "SpatialGridDataFrame") , silent = TRUE)
        if(!class(.Last.value)[1]=="try-error"){
          sel <- m@data[,1]==1 & ( !m@data[,3]==100 | (!m@data[,2]==0|is.na(m@data[,2])) )
          if(sum(sel, na.rm=TRUE)>0){
            m$d <- ifelse(sel, 1, NA) ## | is.na(m@data[,3])
            writeGDAL(m["d"], nfile, type="Byte", mvFlag=0, options="COMPRESS=DEFLATE")
          }
        }
      }             
    }
  }
}
sfInit(parallel=TRUE, cpus=7)
sfLibrary(gdalUtils)
sfLibrary(sp)
sfExport("equi7t3", "tile2.tif", "gdalwarp")
x <- sfLapply(equi7t3, function(t){try(tile2.tif(t))})
sfStop() ## end of parallel processing
tile2.tif(t=equi7t3[[4]])

## Mosaics (per continent):
for(j in 1:length(equi7t3)){
  if(!file.exists(paste0("SMK_", names(equi7t3)[j], "_250m.tif"))){
    t.lst <- list.files(path="stiled", pattern=glob2rx(paste0("SMK_", names(equi7t3)[j], "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    unlink("my_liste.txt")
    cat(t.lst, sep="\n", file="my_liste.txt")
    gdalbuildvrt(input_file_list="my_liste.txt", output.vrt=paste0(names(equi7t3)[j], ".vrt"))
    system(paste0(gdalwarp, ' ', paste0(names(equi7t3)[j], ".vrt"), ' ', paste0("SMK_", names(equi7t3)[j], "_250m.tif"), ' -r \"near\" -srcnodata 0 -dstnodata 0 -ot \"Byte\"'))
  }
}

## LAND MASK
for(j in 1:length(equi7t3)){
  if(!file.exists(paste0("LMK_", names(equi7t3)[j], "_250m.tif"))){
    t.lst <- list.files(path="stiled", pattern=glob2rx(paste0("LMK_", names(equi7t3)[j], "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    unlink("my_liste.txt")
    cat(t.lst, sep="\n", file="my_liste.txt")
    gdalbuildvrt(input_file_list="my_liste.txt", output.vrt=paste0(names(equi7t3)[j], ".vrt"))
    system(paste0(gdalwarp, ' ', paste0(names(equi7t3)[j], ".vrt"), ' ', paste0("LMK_", names(equi7t3)[j], "_250m.tif"), ' -r \"near\" -srcnodata 0 -dstnodata 0 -ot \"Byte\"'))
    gzip(paste0("LMK_", names(equi7t3)[j], "_250m.tif"), overwrite=TRUE)
  }
}
## 1km:
unlink("LMKGLC3a.tif")
system(paste0(gdalwarp, ' lmask.vrt LMKGLC3a.tif -t_srs \"+proj=longlat +datum=WGS84\" -r \"near\" -tr 0.008333333 0.008333333 -te -180 -90 180 90 -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))
unlink("LMKGLC4a.tif")
system(paste0(gdalwarp, ' lmask.vrt LMKGLC4a.tif -t_srs \"+proj=longlat +datum=WGS84\" -r \"near\" -tr 0.004166667 0.004166667 -te -180 -90 180 90 -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))

## GLC classes:
tile_glc.tif <- function(t, t_srs = proj4string(t), cl){
  for(j in 1:nrow(t)){
    nfile <- paste0("tiled/GLC", cl, "_", strsplit(paste(t@data[j,"SHORTNAME"]), " ")[[1]][2], "_", t@data[j,"TILE"], ".tif")
    if(!file.exists(nfile)){
      te <- as.vector(bbox(t[j,]))
      system(paste0(gdalwarp, ' ', paste0('glc', cl, '.vrt'), ' ', nfile, ' -t_srs \"', t_srs, '\" -tr 250 250 -srcnodata 255 -dstnodata 255 -te ', paste(te, collapse=" "), ' -ot \"Byte\" -co \"COMPRESS=DEFLATE\"'))
    }
  }
}
## located tiles with all missing pixels?
#tif.lst <- list.files(path="./tiled", pattern=glob2rx("GLC??_??_???_???.tif$"), full.names=TRUE)

for(i in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
  if(!file.exists(paste0("glc", i, ".vrt"))){
    glc00.lst <- list.files(path="tmp/", pattern=glob2rx(paste0("*_", i,"_ll.tif$")), full.names=TRUE)
    unlink("my_liste.txt")
    cat(glc00.lst, sep="\n", file="my_liste.txt")
    gdalbuildvrt(input_file_list="my_liste.txt", output.vrt=paste0("glc", i, ".vrt"))
  }
  sfInit(parallel=TRUE, cpus=9)
  sfLibrary(gdalUtils)
  sfLibrary(sp)
  sfExport("equi7t3", "tile_glc.tif", "gdalwarp", "i")
  x <- sfLapply(equi7t3, function(x){ tile_glc.tif(t=x, cl=i) })
  sfStop()
}

## Mosaics per continent:
for(cl in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)){
  for(j in 1:length(equi7t3)){
    if(!file.exists(paste0("GLC", cl, "_", names(equi7t3)[j], "_250m.tif.gz"))){
      t.lst <- list.files(path="tiled", pattern=glob2rx(paste0("GLC", cl, "_", names(equi7t3)[j], "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
      unlink("my_liste.txt")
      cat(t.lst, sep="\n", file="my_liste.txt")
      gdalbuildvrt(input_file_list="my_liste.txt", output.vrt=paste0("glc_", names(equi7t3)[j], ".vrt"))
      system(paste0(gdalwarp, ' ', paste0("glc_", names(equi7t3)[j], ".vrt"), ' ', paste0("GLC", cl, "_", names(equi7t3)[j], "_250m.tif"), ' -r \"near\" -srcnodata 255 -dstnodata 255 -ot \"Byte\"'))
      gzip(paste0("GLC", cl, "_", names(equi7t3)[j], "_250m.tif"))
    }
  }
}

## In Sinusoidal projection:
load("E:\\WORLDGRIDS\\MODIS_1km\\modis_sinusoidal\\modis_grid.rda")

tile_sin.tif <- function(t, t_srs = "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs", cl){
  nfile <- paste0("tiled/GLC", cl, "_", modis_grid$TILE[t], "_sin.tif")
  if(!file.exists(nfile)){
    te <- as.vector(bbox(modis_grid[t,]))
    system(paste0(gdalwarp, ' ', paste0('glc', cl, '.vrt'), ' ', nfile, ' -t_srs \"', t_srs, '\" -tr 250 250 -te ', paste(te, collapse=" "), ' -ot \"Byte\"'))
  }
}

## 1 km resolution:
for(i in 1:length(classes.t)){
  nfile <- paste0("G", classes.t[i], "NSD3a.tif")
  system(paste0(gdalwarp, ' ', paste0('glc', classes[i], '.vrt'), ' ', nfile,' -t_srs \"+proj=longlat +datum=WGS84\" -r \"average\" -tr 0.008333333 0.008333333 -te -180 -90 180 90 -ot \"Byte\"'))
}

## per class:
sfInit(parallel=TRUE, cpus=10)
sfLibrary(gdalUtils)
sfLibrary(sp)
sfExport("modis_grid", "tile_sin.tif", "gdalwarp")
x <- sfLapply(1:nrow(modis_grid), function(x){ tile_sin.tif(t=x, cl=50) })
sfStop()

GLC_50.sin <- list.files(path="./tiled", pattern=glob2rx("GLC50_*_sin.tif$"), full.names=TRUE)
unlink("my_liste.txt")
cat(GLC_50.sin, sep="\n", file="my_liste.txt")
gdalbuildvrt(input_file_list="my_liste.txt", output.vrt="glc_sin50.vrt")
system(paste0(gdalwarp, ' glc_sin50.vrt GLC50_sin250m.tif -tr 250 250 -srcnodata 255 -dstnodata 255 -r \"near\" -ot \"Byte\"'))

## Check the boundary line:
ogr2ogr("E:\\WORLDGRIDS\\NaturalEarth\\ne_10m_coastline.shp", "EU_ne_10m_coastline.shp", t_srs=proj4string(equi7t3[["EU"]]))


## 100 m resolution images per continent:
GLC_30m <- list.files(pattern=glob2rx("*_30m.tif$"), full.names=TRUE, recursive = TRUE)
str(GLC_30m) ## 852 tile
## Resample to 100 m:
sfInit(parallel=TRUE, cpus=6)
sfExport("GLC_30m", "gdalwarp")
t <- sfLapply(GLC_30m, function(x){ if(!file.exists(paste0("tmp/", gsub("_30m.tif", "_ll.tif", basename(x))))){ system(paste0(gdalwarp, ' ', x, ' -t_srs \"+proj=longlat +datum=WGS84\" ', paste0("tmp/", gsub("_30m.tif", "_ll.tif", basename(x))), ' -tr 0.0008333333 0.0008333333 -r \"near\" -ot \"Byte\" -dstnodata 255 -co \"COMPRESS=DEFLATE\"')) } })
sfStop()

GLC_30m_ll <- paste0("./tmp/", gsub("_30m.tif", "_ll.tif", basename(GLC_30m)))
unlink("my_liste.txt")
cat(GLC_30m_ll, sep="\n", file="my_liste.txt")
gdalbuildvrt(input_file_list="my_liste.txt", output.vrt="GLC_30m.vrt")

for(j in 1:length(equi7t3)){
  if(!file.exists(paste0("GLC", "_", names(equi7t3)[j], "_100m.tif"))){
    t_srs <- proj4string(equi7t3[[j]])
    te <- as.vector(bbox(equi7t3[[j]]))
    system(paste0(gdalwarp, ' GLC_30m.vrt ', paste0("GLC2010", "_", names(equi7t3)[j], "_100m.tif"), ' -r \"near\" -te ', paste(te, collapse=" "), ' -t_srs \"', t_srs, '\" -tr 100 100 -srcnodata 255 -dstnodata 255 -ot \"Byte\"')) ## 
  }
}

##gzip(paste0("GLC", "_", names(equi7t3)[j], "_30m.tif"))

## end of script;