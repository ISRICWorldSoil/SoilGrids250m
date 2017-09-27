## Functions to calculate uncertainty/errors in SoilGrids250m
## Tom.Hengl@isric.org

setwd("/data/models")
list.of.packages = c("entropy", "plyr", "rgdal")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(entropy)
library(plyr)
library(rgdal)
library(snowfall)
library(raster)
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate =  "gdal_translate"
  gdalwarp =  "gdalwarp"
  gdalbuildvrt = "gdalbuildvrt"
  saga_cmd = "/usr/local/bin/saga_cmd"
}

## Factor type variables
## derive Scaled Shannon Entropy (100 is a maximum error; 0 is perfect prediction)
entropy_tile <- function(i, in.path, varn, levs){
  out.p <- paste0(in.path, "/", i, "/SSI_", varn, "_", i, ".tif")
  if(!file.exists(out.p)){
    tif.lst <- paste0(in.path, "/", i, "/", varn, "_", levs, "_", i, ".tif")
    s <- raster::stack(tif.lst)
    s <- as(as(s, "SpatialGridDataFrame"), "SpatialPixelsDataFrame")
    gc()
    v <- unlist(alply(s@data, 1, .fun=function(x){entropy.empirical(unlist(x))})) 
    s$SSI <- round(v/entropy.empirical(rep(1/length(levs),length(levs)))*100)
    writeGDAL(s["SSI"], out.p, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
  }
}

## Run in parallel (VERY COMPUTATIONALY INTENSIVE):
pr.dirs <- basename(list.dirs("/data/tt/SoilGrids250m/predicted250m")[-1])

varn = "TAXNWRB"
levs = list.files(path="/data/tt/SoilGrids250m/predicted250m/AF_006_070", pattern=glob2rx(paste0(varn, "_*_AF_*_*.tif$")))
levs = sapply(levs, function(x){strsplit(x, "_")[[1]][2]})
sfInit(parallel=TRUE, cpus=48)
sfExport("entropy_tile", "levs", "varn")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(entropy)
out <- sfClusterApplyLB(pr.dirs, function(i){try( entropy_tile(i, in.path="/data/tt/SoilGrids250m/predicted250m", varn, levs) )})
sfStop()

varn = "TAXOUSDA"
levs = list.files(path="/data/tt/SoilGrids250m/predicted250m/AF_006_070", pattern=glob2rx(paste0(varn, "_*_AF_*_*.tif$")))
levs = sapply(levs, function(x){strsplit(x, "_")[[1]][2]})
sfInit(parallel=TRUE, cpus=48)
sfExport("entropy_tile", "levs", "varn")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(entropy)
out <- sfClusterApplyLB(pr.dirs, function(i){try( entropy_tile(i, in.path="/data/tt/SoilGrids250m/predicted250m", varn, levs) )})
sfStop()

## clean up:
#del.lst = list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("SSI*.tif"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## 1 km resolution
obj <- GDALinfo("/data/GEOG/TAXOUSDA_1km_ll.tif")
r1km = raster("/data/GEOG/TAXOUSDA_1km_ll.tif")
tile.lst <- getSpatialTiles(obj, block.x=5, return.SpatialPolygons=TRUE)
tile.tbl <- getSpatialTiles(obj, block.x=5, return.SpatialPolygons=FALSE)
tile.tbl$ID = as.character(1:nrow(tile.tbl))
## 2160 tiles
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
writeOGR(tile.pol, "tiles_5deg.shp", "tiles_5deg", "ESRI Shapefile")
system(paste('gdal_translate /data/GEOG/TAXOUSDA_1km_ll.tif /data/tiled1km/TAXOUSDA_1km_ll.sdat -of \"SAGA\" -ot \"Byte\"'))
system(paste0(saga_cmd, ' -c=48 shapes_grid 2 -GRIDS=\"/data/tiled1km/TAXOUSDA_1km_ll.sgrd\" -POLYGONS=\"tiles_5deg.shp\" -PARALLELIZED=1 -RESULT=\"ov_ADMIN_tiles_5deg.shp\"'))
ov_ADMIN_5deg = readOGR("ov_ADMIN_tiles_5deg.shp", "ov_ADMIN_tiles_5deg")
summary(sel.blocks <- !is.na(ov_ADMIN_5deg$TAXOUSDA_1k.5))

entropy_block <- function(i, tile.tbl, tif.lst){
  out.p <- paste0("/data/tiled1km/SSI_", strsplit(basename(tif.lst), "_")[[1]][1], "_", i, ".tif")
  if(!file.exists(out.p)){
    s = readGDAL(fname=tif.lst[1], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    s <- as(s, "SpatialPixelsDataFrame")
    for(j in 2:length(tif.lst)){
      s@data[,j] = readGDAL(fname=tif.lst[j], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[s@grid.index]
    }
    gc()
    v <- unlist(alply(s@data, 1, .fun=function(x){entropy.empirical(unlist(x))})) 
    s$SSI <- round(v/entropy.empirical(rep(1/length(tif.lst),length(tif.lst)))*100)
    #plot(raster(s["SSI"]), col=SAGA_pal[[1]])
    writeGDAL(s["SSI"], out.p, type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
  }
}

## TAKES CA 1 hr
varn = "TAXOUSDA"
levs = list.files(path="/data/tt/SoilGrids250m/predicted250m/AF_006_070", pattern=glob2rx(paste0(varn, "_*_AF_*_*.tif$")))
levs = sapply(levs, function(x){strsplit(x, "_")[[1]][2]})
tif.lst <- paste0("/data/GEOG/", varn, "_", levs, "_1km_ll.tif")
any(!file.exists(tif.lst))
sfInit(parallel=TRUE, cpus=48)
sfExport("entropy_block", "tile.tbl", "tif.lst")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(entropy)
out <- sfClusterApplyLB(which(sel.blocks), function(i){try( entropy_block(i, tile.tbl, tif.lst) )})
sfStop()

## TAKES CA 2-3 hrs because number of classes is much larger!
varn2 = "TAXNWRB"
levs2 = list.files(path="/data/tt/SoilGrids250m/predicted250m/AF_006_070", pattern=glob2rx(paste0(varn2, "_*_AF_*_*.tif$")))
levs2 = sapply(levs2, function(x){strsplit(x, "_")[[1]][2]})
tif.lst2 <- paste0("/data/GEOG/", varn2, "_", levs2, "_1km_ll.tif")
any(!file.exists(tif.lst2))
sfInit(parallel=TRUE, cpus=48)
sfExport("entropy_block", "tile.tbl", "tif.lst2")
sfLibrary(raster)
sfLibrary(rgdal)
sfLibrary(plyr)
sfLibrary(entropy)
out <- sfClusterApplyLB(which(sel.blocks), function(i){try( entropy_block(i, tile.tbl, tif.lst2) )})
sfStop()

tmp.lst <- list.files(path="/data/tiled1km", pattern=glob2rx(paste0("SSI_TAXNWRB_*.tif$")), full.names=TRUE, recursive=TRUE)
out.tmp <- tempfile(fileext = ".txt")
vrt.tmp <- tempfile(fileext = ".vrt")
cat(tmp.lst, sep="\n", file=out.tmp)
system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
system(paste0(gdalwarp, ' ', vrt.tmp, ' /data/GEOG/SSE_TAXNWRB_1km_ll.tif -ot \"Byte\" -dstnodata \"255\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -r \"near\"'))

tmp.lst <- list.files(path="/data/tiled1km", pattern=glob2rx(paste0("SSI_TAXOUSDA_*.tif$")), full.names=TRUE, recursive=TRUE)
out.tmp <- tempfile(fileext = ".txt")
vrt.tmp <- tempfile(fileext = ".vrt")
cat(tmp.lst, sep="\n", file=out.tmp)
system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
system(paste0(gdalwarp, ' ', vrt.tmp, ' /data/GEOG/SSE_TAXOUSDA_1km_ll.tif -ot \"Byte\" -dstnodata \"255\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -r \"near\"'))
