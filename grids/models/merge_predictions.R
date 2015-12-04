## Create mosaicks from predictions (ca. 2500 EQUI7T3 tiles) - SoilGrids250m
## Tom.Hengl@isric.org

library(snowfall)
library(rgdal)
gdalwarp = "/usr/local/bin/gdalwarp"
gdalbuildvrt = "/usr/local/bin/gdalbuildvrt"
gdal_translate =  "/usr/local/bin/gdal_translate"
gdal_merge.py = "/usr/local/bin/gdal_merge.py"
system("/usr/local/bin/gdal-config --version")
load("equi7t3.rda")

## Create dirs:
#x <- lapply(paste0("/data/GEOG/", names(equi7t3)), dir.create, recursive=TRUE)

## clean-up:
#del.lst <- list.files(path="/data/GEOG", pattern=glob2rx("*_ll.tif$"), full.names=TRUE, recursive=TRUE)
#unlink(del.lst)

## Create mosaicks:
mosaick.equi7t3 <- function(i, j, varn, in.path, r, te, tr){
  if(i=="dominant"){
    out.tif <- paste0('/data/GEOG/', j, '/', varn, '_', j, '_250m_ll.tif')
  } else {
    out.tif <- paste0('/data/GEOG/', j, '/', varn, '_', i, '_', j, '_250m_ll.tif')
  }
  if(!file.exists(out.tif)){
    if(i=="dominant"){
      tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    } else {
      tmp.lst <- list.files(path=in.path, pattern=glob2rx(paste0(varn, "_", i, "_", j, "_*_*.tif$")), full.names=TRUE, recursive=TRUE)
    }
    if(length(tmp.lst)>0){
      out.tmp <- tempfile(fileext = ".txt")
      vrt.tmp <- tempfile(fileext = ".vrt")
      cat(tmp.lst, sep="\n", file=out.tmp)
      system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
      system(paste0(gdalwarp, ' ', vrt.tmp, ' ', out.tif, ' -t_srs \"+proj=longlat +datum=WGS84\" -r \"', r,'\" -ot \"Byte\" -dstnodata \"255\" -te ', paste(te, collapse=" "),' -tr ', tr, ' ', tr, ' -co \"BIGTIFF=YES\" -wm 4000')) ## -co \"COMPRESS=DEFLATE\"
    }
  }
}

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

## test it
#x = sapply(1:length(equi7t3), function(x){mosaick.equi7t3(j=names(equi7t3)[x], i="Fibrists", varn="TAXOUSDA", te=ext[[x]])})

## Reprojection TAKES CA 3-5 hrs (compression is the most time-consuming?)
#m <- readRDS("/data/models/TAXOUSDA/m_TAXOUSDA.rds")
#levs <- m_TAXOUSDA$lev
#m <- readRDS("/data/models/TAXNWRB/m_TAXNWRB.rds")
#levs <- gsub(" ", "\\.", gsub("\\)", "\\.", gsub(" \\(", "\\.\\.", m$lev)))
#rm(m)

## Merge everything into a single mosaick
make_mosaick <- function(i, varn, ext, resample1="bilinear", resample2="average", r250m = TRUE, tr = 0.002083333, r="bilinear", in.path="/data/predicted"){
  if(i=="dominant"){
    out.tif <- paste0("/data/GEOG/", varn, "_250m_ll.tif")
  } else {
    out.tif <- paste0("/data/GEOG/", varn, "_", i, "_250m_ll.tif")
  }
  if(!file.exists(out.tif)){
    ## per continent:
    x <- sapply(1:length(equi7t3), function(x){mosaick.equi7t3(j=names(equi7t3)[x], i=i, varn=varn, te=ext[[x]], tr=tr, r=r, in.path=in.path)})
    if(i=="dominant"){
      in.tif <- list.files(path="/data/GEOG", pattern=glob2rx(paste0(varn, "_*_250m_ll.tif$")), full.names=TRUE, recursive=TRUE)
    } else {
      in.tif <- list.files(path="/data/GEOG", pattern=glob2rx(paste0(varn, "_", i, "_*_250m_ll.tif$")), full.names=TRUE, recursive=TRUE)
    }
    if(length(in.tif)>0){
      out.tmp <- tempfile(fileext = ".txt")
      vrt.tmp <- tempfile(fileext = ".vrt")
      cat(in.tif[c(2,3,5,7,6,1,4)], sep="\n", file=out.tmp)
      system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
      ## relatively fast
      if(r250m == TRUE){
        system(paste0(gdal_translate, ' -of GTiff ', vrt.tmp, ' ', out.tif, ' -ot \"Byte\" -a_nodata \"255\" -r \"', resample1, '\" -co \"COMPRESS=DEFLATE\" -co \"TILED=YES\" -co \"BLOCKXSIZE=512\" -co \"BLOCKYSIZE=512\" -co \"BIGTIFF=YES\"'))
      }
      system(paste0(gdal_translate, ' -of GTiff -r \"', resample2,'\" -tr 0.008333333 0.008333333 ', vrt.tmp, ' ', gsub("250m_ll.tif", "1km_ll.tif", out.tif), ' -ot \"Byte\" -a_nodata \"255\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\"'))
      unlink(out.tmp)
      unlink(vrt.tmp)
      unlink(in.tif)      
    }
  }
}

## Parellel mosaicks per continent:
sfInit(parallel=TRUE, cpus=40)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdal_translate", "ext", "levs", "mosaick.equi7t3", "make_mosaick")
#out <- sfClusterApplyLB(levs, function(i){make_mosaick(i, varn="TAXOUSDA", ext=ext)})
out <- sfClusterApplyLB(levs, function(i){make_mosaick(i, varn="TAXNWRB", ext=ext)})
sfStop()

## only dominant class:
make_mosaick(i="dominant", varn="TAXNWRB", ext=ext, resample1="near", resample2="near", r="near")
make_mosaick(i="dominant", varn="TAXOUSDA", ext=ext, resample1="near", resample2="near", r="near")

## test soil property:
#make_mosaick(i="M_sd1", varn="ORCDRC", ext=ext, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE)

## Soil depths:
tbdr = c("BDRICM", "BDRLOG", "BDTICM")
sfInit(parallel=TRUE, cpus=3)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdal_translate", "ext", "tbdr", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(tbdr, function(x){make_mosaick(i="M", varn=x, ext=ext, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE)})
sfStop()

## Resample all soil properties:
tvars = c("ORCDRC", "PHIHOX", "CRFVOL", "SNDPPT", "SLTPPT", "CLYPPT", "BLD", "CECSUM")
props = rep(tvars, 6)
varn.lst =  paste0("M_sd", sapply(1:6, function(x){rep(x, length(tvars))}))
sfInit(parallel=TRUE, cpus=40)
sfExport("equi7t3", "gdalbuildvrt", "gdalwarp", "gdal_translate", "ext", "props", "varn.lst", "mosaick.equi7t3", "make_mosaick")
out <- sfClusterApplyLB(1:length(props), function(x){make_mosaick(varn.lst[x], varn=props[x], ext=ext, tr=0.008333333, in.path="/data/predicted1km", r250m=FALSE)})
sfStop()
