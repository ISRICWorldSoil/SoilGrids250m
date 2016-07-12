## Create mosaicks from predictions (ca. 2500 EQUI7T3 tiles) - SoilGrids250m
## Tom.Hengl@isric.org

setwd("/data/models")
library(snowfall)
library(rgdal)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
gdal_translate =  "/usr/local/bin/gdal_translate"
gdal_merge.py = "/usr/local/bin/gdal_merge.py"
gdaladdo = "/usr/local/bin/gdaladdo"
system("/usr/local/bin/gdal-config --version")
load("equi7t3.rda")
load("equi7t1.rda")
source("mosaick_functions.R")
des <- read.csv("SoilGrids250m_COVS250m.csv")
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

## Reprojection TAKES CA 3-5 hrs
tvars = c("ORCDRC", "PHIHOX", "PHIKCL", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "CECSOL", "BLDFIE", "TEXMHT", "AWCh1", "AWCh2", "AWCh3", "WWP", "AWCtS")
#tvars = c("ORCDRC", "PHIHOX", "PHIKCL", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "CECSUM", "BLD", "TEXMHT")
props = c(rep(tvars, 7), rep("OCSTHA", 6))
varn.lst = c(paste0("M_sl", sapply(1:7, function(x){rep(x, length(tvars))})), paste0("M_sd", 1:6))
ot.lst <- c(rep(c("Int16","Byte","Byte","Byte","Byte","Byte","Byte","Int16","Int16","Byte","Byte","Byte","Byte","Byte","Byte"), 7), rep("Int16", 6))
dstnodata.lst <- c(rep(c(-32768, 255, 255, 255, 255, 255, 255, -32768, -32768, 255, 255, 255, 255, 255, 255), 7), rep(-32768, 6))
resample1.lst <- c(rep(c("bilinear", "bilinear", "bilinear", "bilinear", "bilinear", "bilinear", "bilinear", "bilinear", "bilinear", "near", "bilinear", "bilinear", "bilinear", "bilinear", "bilinear"), 7), rep("bilinear", 6))
resample2.lst <- c(rep(c("average", "average", "average", "average", "average", "average", "average", "average", "average", "mode", "average", "average", "average", "average", "average"), 7), rep("average", 6))

## test soil property:
#make_mosaick(i="M_sd1", varn="ORCDRC", ext=ext, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE)

## Resample all soil properties to 250m:
## TAKES >7-8 hours
## without compression TAKES A LOT OF HARD DISK SPACE
sfInit(parallel=TRUE, cpus=28)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "tile.names", "gdaladdo", "gdal_translate", "ext", "props", "varn.lst", "mosaick.equi7t3", "make_mosaick", "ot.lst", "dstnodata.lst", "resample1.lst", "resample2.lst")
out <- sfClusterApplyLB(1:length(props), function(x){make_mosaick(varn.lst[x], varn=props[x], ext=ext, in.path="/data/predicted", ot=ot.lst[x], dstnodata=dstnodata.lst[x], tile.names=tile.names, tr=0.002083333, r250m=TRUE, r=resample1.lst[x], resample1=resample1.lst[x], resample2=resample2.lst[x])})
sfStop()

## Rename some soil properties:
from.lst = c("CECSUM_M","BLD_M")
to.lst = c("CECSOL_M", "BLDFIE_M")
for(j in 1:length(from.lst)){
  from <- list.files(path="/data/GEOG", pattern=glob2rx(paste0("^",from.lst[j],"*.tif")), full.names = TRUE)
  to <- gsub(from.lst[j], to.lst[j], from) 
  x= sapply(1:length(from), function(x){ if(file.exists(from[x])){ file.rename(from=from[x], to=to[x]) } } )
}

## Parellel mosaicks per continent:
m <- readRDS("/data/models/TAXNWRB/mnetX_TAXNWRB.rds")
levs <- gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", m$finalModel$lev)))
sfInit(parallel=TRUE, cpus=48)
sfExport("gdalbuildvrt", "gdalwarp", "gdaladdo", "tile.names", "gdal_translate", "ext", "levs", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(levs, function(i){make_mosaick(i, varn="TAXNWRB", ext=ext, tile.names=tile.names, in.path="/data/predicted", tr=0.002083333, r250m=TRUE)})
#out <- sfClusterApplyLB(levs, function(i){make_mosaick(i, varn="TAXNWRB", ext=ext, tile.names=tile.names, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE)})
sfStop()

m <- readRDS("/data/models/TAXOUSDA/mnetX_TAXOUSDA.rds")
levs <- m$finalModel$lev
rm(m)
sfInit(parallel=TRUE, cpus=48)
sfExport("gdalbuildvrt", "gdalwarp", "gdaladdo", "tile.names", "gdal_translate", "ext", "levs", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(levs, function(i){make_mosaick(i, varn="TAXOUSDA", ext=ext, tile.names=tile.names, in.path="/data/predicted", tr=0.002083333, r250m=TRUE)})
#out <- sfClusterApplyLB(levs, function(i){make_mosaick(i, varn="TAXOUSDA", ext=ext, tile.names=tile.names, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE)})
sfStop()

## Land mask:
#make_mosaick(i="dominant", varn="LMK", ext=ext, resample1="near", resample2="near", r="near", in.path="/data/mask", tile.names=tile.names, tr=0.002083333, r250m=TRUE, ot="Byte")
#make_mosaick(i="dominant", varn="TAXOUSDA", ext=ext, resample1="near", resample2="mode", r="near", tr=0.002083333, r250m=TRUE, tile.names=tile.names)

## only dominant class:
varn.cl.lst <- c("TAXNWRB","TAXOUSDA")
sfInit(parallel=TRUE, cpus=length(varn.cl.lst))
sfExport("gdalbuildvrt", "gdalwarp", "gdaladdo", "tile.names", "gdal_translate", "ext",  "mosaick.equi7t3", "varn.cl.lst", "make_mosaick")
out <- sfClusterApplyLB(1:length(varn.cl.lst), function(x){ make_mosaick(i="dominant", varn=varn.cl.lst[x], ext=ext, resample1="near", resample2="mode", r="near", tr=0.002083333, r250m=TRUE, tile.names=tile.names) })
sfStop()

## Organic soils:
make_mosaick(i="dominant", varn="HISTPR", ext=ext, in.path="/data/predicted", tr=0.002083333, r250m=TRUE, ot="Int16", dstnodata=-32768, tile.names=tile.names)


## 1km
sfInit(parallel=TRUE, cpus=48)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "ext", "props", "varn.lst", "mosaick.equi7t3", "make_mosaick", "tile.names")
out <- sfClusterApplyLB(1:length(props), function(x){make_mosaick(varn.lst[x], varn=props[x], ext=ext, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE, ot="Int16", dstnodata=-32768, tile.names=tile.names)})
sfStop()

tbdr = c("BDRICM", "BDRLOG", "BDTICM")
## Soil depths 250m:
#make_mosaick(i="M", varn="BDTICM", ext=ext, tr=0.002083333, in.path="/data/predicted", r250m=TRUE, ot="Int32", dstnodata=-99999, tile.names=tile.names)
sfInit(parallel=TRUE, cpus=3)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "ext", "tbdr", "mosaick.equi7t3", "make_mosaick", "tile.names")
out <- sfClusterApplyLB(tbdr, function(x){try( make_mosaick(i="M", varn=x, ext=ext, in.path="/data/predicted", tr=0.002083333, r250m=TRUE, ot="Int32", dstnodata=-99999, tile.names=tile.names) )})
sfStop()

## Soil depths 1km:
sfInit(parallel=TRUE, cpus=3)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "ext", "tbdr", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(tbdr, function(x){make_mosaick(i="M", varn=x, ext=ext, tr=0.008333333, in.path="/data/predicted", r250m=FALSE, ot="Int32", dstnodata=-99999, tile.names=tile.names)})
sfStop()

## Covariates at 250m:
tcovs <- as.character(des$WORLDGRIDS_CODE[c(grep(glob2rx("***MOD5"), des$WORLDGRIDS_CODE), grep(glob2rx("***MOD4"), des$WORLDGRIDS_CODE))])
sfInit(parallel=TRUE, cpus=ifelse(length(tcovs)>10,10,length(tcovs)))
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdaladdo", "gdal_translate", "ext", "tcovs", "mosaick.equi7t3", "make_mosaick", "tile.names")
out <- sfClusterApplyLB(tcovs, function(x){try( make_mosaick(i="dominant", varn=x, ext=ext, in.path="/data/covs1t", tr=0.002083333, r250m=TRUE, ot="Int16", dstnodata=-32768, tile.names=tile.names) )})
sfStop()
