## Generation of global DEM parameters at 250m resolution
## tom.hengl@isric.org

setwd("/data/MDEM")
load(".RData")
library(R.utils)
library(rgdal)
library(snowfall)
library(pkgmaker)
library(raster)
if(.Platform$OS.type == "windows"){
  gdal.dir <- shortPathName("C:/Program files/GDAL")
  gdal_translate <- paste0(gdal.dir, "/gdal_translate.exe")
  gdalwarp <- paste0(gdal.dir, "/gdalwarp.exe") 
} else {
  gdal_translate = "gdal_translate"
  gdalwarp = "gdalwarp"
  saga_cmd = "saga_cmd"
}
## DEM derivation functions:
source("DEM_functions.R")
source("/data/models/mosaick_functions.R")
## regions of interest:
regs <- c("AF","AS","EU","NA","OC","SA")

## Derive DEM parameters for all continents:
#sapply(regs, function(x){ saga_DEM_derivatives(INPUT=paste0("MDEM_", x, "_250m.tif"), sel=c("CRV","DVM","MRN","TPI")) })
sapply(regs, function(x){ saga_DEM_derivatives(INPUT=paste0("MDEM_", x, "_250m.sgrd"), sel=c("CRV","DVM","MRN","TPI")) })
unlink(paste0("MDEM_", regs, "_250m.sdat"))

## Convert numbers to integers + convert to Geotiffs:
des <- read.csv("/data/models/SoilGrids250m_COVS250m.csv")
in.lst = list.files("/data/MDEM", pattern=glob2rx("MDEM_*_250m_*.sdat$"), full.names = TRUE)

multiply10 = function(x,m){ x*10 }
multiply100 = function(x,m){ x*100 }
multiply1000 = function(x,m){ x*1000 }
## run in loop (takes ca 5 hours)
for(i in in.lst){
  filename = gsub(".sdat", ".tif", i)
  if(!file.exists(filename)){
    s = raster(i)
    beginCluster()
    if(length(grep("devmean",i))>0|length(grep("tpi",i))>0|length(grep("slope",i))>0|length(grep("vbf",i))>0){
      r0 <- clusterR(s, fun=calc, args=list(fun=multiply100), filename=filename, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
    }
    if(length(grep("down",i))>0|length(grep("uplocal",i))>0|length(grep("open",i))>0){
      r0 <- clusterR(s, fun=calc, args=list(fun=multiply1000), filename=filename, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
    }
    if(length(grep("mrn",i))>0|length(grep("TWI",i))>0|length(grep("vdepth",i))>0){
      r0 <- clusterR(s, fun=calc, args=list(fun=multiply10), filename=filename, datatype="INT2S", options=c("COMPRESS=DEFLATE"))
    }
    endCluster()
  }
}

## clean-up
unlink(in.lst)

## grid definition:
r <- raster("/data/GEOG/TAXOUSDA_250m_ll.tif")
ncols = ncol(r)
nrows = nrow(r)
cellsize = res(r)[1]

## Derive global mosaics (DEM derivatives)
dem.lst = c("devmean","devmean2","tpi","mrn","twi","downlocal","down","uplocal")
out.dem.lst = c("DVMMRG5", "DV2MRG5", "TPIMRG5", "MRNMRG5", "TWIMRG5", "CRDMRG5", "CRVMRG5", "CRUMRG5") #, "NEGMRG5", "POSMRG5", "SLPMRG5", , "VBFMRG5")
tiled.lst = lapply(dem.lst, function(x){list.files(path="/data/MDEM", pattern=glob2rx(paste0("MDEM_*_250m_", x, ".tif$")), full.names=TRUE)})
names(tiled.lst) = out.dem.lst
system(paste0('gdalinfo ', tiled.lst[[1]][2])) 
nodata.DEM <- des$NO_DATA[match(out.dem.lst, des$WORLDGRIDS_CODE)]
ot.DEM <- as.character(des$DATA_FORMAT[match(out.dem.lst, des$WORLDGRIDS_CODE)])
tile.names = sapply(basename(tiled.lst[[x]]), function(i){strsplit(i, "_")[[1]][2]})
ext.DEM <- ext[tile.names]
## in parallel:
sfInit(parallel=TRUE, cpus=length(out.dem.lst))
sfExport("equi7t1", "ext.DEM", "out.dem.lst", "tiled.lst", "mosaick.equi7t3", "make_mosaick", "nodata.DEM", "ot.DEM", "tile.names")
out <- sfClusterApplyLB(1:length(out.dem.lst), function(x){try( make_mosaick(i="dominant", varn=out.dem.lst[x], ext=ext.DEM, tr=0.002083333, r250m=TRUE, ot=ot.DEM[x], dstnodata=nodata.DEM[x], tile.names=tile.names, build.pyramids=FALSE, vrt.tmp=tiled.lst[[x]], cleanup=FALSE) )})
sfStop()
