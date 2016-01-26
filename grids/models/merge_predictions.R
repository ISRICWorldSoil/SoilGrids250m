## Create mosaicks from predictions (ca. 2500 EQUI7T3 tiles) - SoilGrids250m
## Tom.Hengl@isric.org

library(snowfall)
library(rgdal)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
gdal_translate =  "/usr/local/bin/gdal_translate"
gdal_merge.py = "/usr/local/bin/gdal_merge.py"
gdaladdo = "/usr/local/bin/gdaladdo"
system("/usr/local/bin/gdal-config --version")
load("equi7t3.rda")
source("mosaick_functions.R")
## Continents:
ext <- as.list(1:7)
#names(ext) <- names(equi7t3)
ext[[1]] <- c(-32.42, -42.90, 64.92, 41.08) ## "AF"
ext[[2]] <- c(-66.4, -56.17, 55.44, -44.30) ## "AN"
ext[[3]] <- c(40.35, -4.67, 180, 87.37) ## "AS"
ext[[4]] <- c(-31.4, 32.2, 60.6, 82.40) ## "EU"
ext[[5]] <- c(-180, -9.71, -10.3, 83.3) ## "NA"
ext[[6]] <- c(92.68, -53.25, 180, 26.38) ## "OC"
ext[[7]] <- c(-122.85, -56.29, -16.04, 20.23) ## "SA"
tile.names <- names(equi7t3)

## Create dirs:
#x <- lapply(paste0("/data/GEOG/", tile.names, dir.create, recursive=TRUE)

## clean-up:
#del.lst <- list.files(path="/data/GEOG", pattern=glob2rx("*_ll.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## test it
#x = sapply(1:length(equi7t3), function(x){mosaick.equi7t3(j=tile.names[x], i="Fibrists", r="bilinear", in.path="/data/predicted", varn="TAXOUSDA", te=ext[[x]], ot="Byte", dstnodata=255, tr=0.002083333)})
#make_mosaick(i="Fibrists", varn="TAXOUSDA", ext=ext, tr=0.002083333, in.path="/data/predicted", r250m=TRUE, tile.names=tile.names)

## Reprojection TAKES CA 3-5 hrs (compression is the most time-consuming?)

## Parellel mosaicks per continent:
m <- readRDS("/data/models/TAXNWRB/m_TAXNWRB.rds")
levs <- gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", m$lev)))
sfInit(parallel=TRUE, cpus=40)
sfExport("gdalbuildvrt", "gdalwarp", "gdaladdo", "tile.names", "gdal_translate", "ext", "levs", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(levs, function(i){make_mosaick(i, varn="TAXNWRB", ext=ext, tile.names=tile.names, in.path="/data/predicted", tr=0.002083333, r250m=TRUE)})
sfStop()

m <- readRDS("/data/models/TAXOUSDA/m_TAXOUSDA.rds")
levs <- m$lev
rm(m)
sfInit(parallel=TRUE, cpus=40)
sfExport("gdalbuildvrt", "gdalwarp", "gdaladdo", "tile.names", "gdal_translate", "ext", "levs", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(levs, function(i){make_mosaick(i, varn="TAXOUSDA", ext=ext, tile.names=tile.names, in.path="/data/predicted", tr=0.002083333, r250m=TRUE)})
sfStop()

## Land mask:
make_mosaick(i="dominant", varn="LMK", ext=ext, resample1="near", resample2="near", r="near", in.path="/data/mask", tile.names=tile.names, tr=0.002083333, r250m=TRUE, ot="Byte")

## only dominant class:
make_mosaick(i="dominant", varn="TAXNWRB", ext=ext, resample1="near", resample2="near", r="near", tile.names=tile.names)
make_mosaick(i="dominant", varn="TAXOUSDA", ext=ext, resample1="near", resample2="near", r="near", tile.names=tile.names)
## Organic soils:
make_mosaick(i="dominant", varn="HISTPR", ext=ext, tile.names=tile.names, ot="Int16", dstnodata=-32768)

## test soil property:
#make_mosaick(i="M_sd1", varn="ORCDRC", ext=ext, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE)

## Soil depths:
tbdr = c("BDRICM", "BDRLOG", "BDTICM")
sfInit(parallel=TRUE, cpus=3)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "ext", "tbdr", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(tbdr, function(x){make_mosaick(i="M", varn=x, ext=ext, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE, ot="Int32", dstnodata=-99999, tile.names=tile.names)})
sfStop()

## Resample all soil properties 1km:
# tvars = c("ORCDRC", "PHIHOX", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "CECSUM", "OCSTHA")
# props = rep(tvars, 6)
# varn.lst =  paste0("M_sd", sapply(1:6, function(x){rep(x, length(tvars))}))
# sfInit(parallel=TRUE, cpus=40)
# sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "ext", "props", "varn.lst", "mosaick.equi7t3", "make_mosaick")
# out <- sfClusterApplyLB(1:length(props), function(x){make_mosaick(varn.lst[x], varn=props[x], ext=ext, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE, ot="Int16", dstnodata=-32768, tile.names=tile.names)})
# sfStop()

## Resample all soil properties to 250m:
## TAKES >7-8 hours
tvars = c("ORCDRC", "PHIHOX", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "CECSUM") ## , "OCSTHA"
props = rep(tvars, 6)
varn.lst =  paste0("M_sd", sapply(1:6, function(x){rep(x, length(tvars))}))
## TAKES A LOT OF HARD DISK SPACE
sfInit(parallel=TRUE, cpus=40)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "tile.names", "gdaladdo", "gdal_translate", "ext", "props", "varn.lst", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(1:length(props), function(x){make_mosaick(varn.lst[x], varn=props[x], ext=ext, in.path="/data/predicted", ot="Int16", dstnodata=-32768, tile.names=tile.names, tr=0.002083333, r250m=TRUE, compress=TRUE)})
sfStop()
