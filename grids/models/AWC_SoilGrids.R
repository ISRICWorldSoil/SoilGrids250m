## Deriving Available Water capacity using SoilGrids250m
## prepared by: T. Hengl (tom.hengl@isric.org)

library(raster)
library(GSIF)
library(rgdal)
library(sp)
library(snowfall)


predictAWC <- function(i, in.path, depths=1:7){
  ## run for each depth:
  for(d in depths){
    filen <- paste0(in.path, "/", i, "/", c("AWCh1", "AWCh2", "AWCh3", "WWP", "AWCtS"), "_M_sl", d, "_", i, ".tif")
    tifn <- c("BLD", "CECSUM", "CLYPPT", "ORCDRC", "PHIHOX", "SLTPPT", "SNDPPT")
    tif.lst <- paste0(in.path, "/", i, "/", tifn, "_M_sl", d, "_", i, ".tif")
    if(any(!file.exists(tif.lst))){ 
      warning(paste("Layers for tile:", i, "missing")) 
    } else {
      if(any(!file.exists(filen))){
        r <- readGDAL(tif.lst[1])
        for(k in 2:length(tif.lst)){
            r@data[,k] <- readGDAL(tif.lst[k])$band1 
        }
        r <- as(r, "SpatialPixelsDataFrame")
        names(r) <- tifn
        r@data <- AWCPTF(r$SNDPPT, r$SLTPPT, r$CLYPPT, r$ORCDRC, r$BLD, r$CECSUM, r$PHIHOX/10)
        for(j in 1:ncol(r@data)){
          r@data[,j] <- round(100*r@data[,j])
          writeGDAL(r[j], filen[j], type="Byte", mvFlag=255, options="COMPRESS=DEFLATE")
        }
      }
    }
  }
}

## test it:
predictAWC(i="NA_060_036", in.path="/data/predicted")

## Run in parallel:
pr.dirs <- basename(list.dirs("/data/predicted")[-1])

sfInit(parallel=TRUE, cpus=48)
sfLibrary(rgdal)
sfLibrary(sp)
sfLibrary(raster)
sfLibrary(GSIF)
sfExport("pr.dirs", "predictAWC")
out <- sfClusterApplyLB(pr.dirs, function(i){try( predictAWC(i, in.path="/data/predicted") )})
sfStop()