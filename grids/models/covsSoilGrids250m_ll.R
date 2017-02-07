## Stacking covariates at 250 m resolution for SoilGrids
## tom.hengl@isric.org

library(RCurl)
library(rgdal)
library(raster)
library(parallel)
library(snowfall)
## Grid definition:
r <- raster("/data/GEOG/TAXOUSDA_250m_ll.tif")
extent(r)
ncols = ncol(r)
nrows = nrow(r)
xllcorner = extent(r)[1]
yllcorner = extent(r)[3]
xurcorner = extent(r)[2]
yurcorner = extent(r)[4]
cellsize = res(r)[1]

des <- read.csv("SoilGrids250m_COVS250m.csv")
## check if all files exist:
des.s = des[!des$ROOT_FILE=="",c("WORLDGRIDS_CODE","ROOT_FILE","RESAMPLE_METHOD")]
s = file.exists(paste(des.s$ROOT_FILE))
summary(s)
#system(paste0('gdalinfo /data/MOD13Q1/M_listeJanFeb.vrt'))

## process in parallel (250m):
sfInit(parallel=TRUE, cpus=ifelse(length(tcovs)>24, 24, length(tcovs)))
sfExport("r", "cellsize", "des.s")
sfLibrary(rgdal)
sfLibrary(raster)
out <- sfClusterApplyLB(1:nrow(des.s), function(k){ if(!file.exists(paste0('/data/stacked250m/', des.s$WORLDGRIDS_CODE[k], '.tif'))){ system(paste0('gdalwarp ', des.s$ROOT_FILE[k], ' /data/stacked250m/', des.s$WORLDGRIDS_CODE[k], '.tif -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -r \"', des.s$RESAMPLE_METHOD[k], '\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "))) } } )
sfStop()

## 1 km
sfInit(parallel=TRUE, cpus=ifelse(length(tcovs)>24, 24, length(tcovs)))
sfExport("r", "des.s")
sfLibrary(rgdal)
sfLibrary(raster)
out <- sfClusterApplyLB(1:nrow(des.s), function(k){ if(!file.exists(paste0('/data/stacked1km/', des.s$WORLDGRIDS_CODE[k], '.tif'))){ system(paste0('gdalwarp ', des.s$ROOT_FILE[k], ' /data/stacked1km/', des.s$WORLDGRIDS_CODE[k], '.tif -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', 1/120, ' ', 1/120, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "))) } } )
sfStop()


## Covariates in Equi7 system (tiles)
library(snowfall)
source("mosaick_functions.R")
sel.equi7 = c(grep(glob2rx("***USG5"), des$WORLDGRIDS_CODE), grep(glob2rx("***MRG5"), des$WORLDGRIDS_CODE))
tcovs <- as.character(des$WORLDGRIDS_CODE[sel.equi7])
nodata.lst <- des$NO_DATA[sel.equi7]
ot.lst <- as.character(des$DATA_FORMAT[sel.equi7])
## test it:
#make_mosaick(i="dominant", varn=tcovs[1], ext=ext, in.path="/data/tt/SoilGrids250m/covs1t", tr=0.002083333, r250m=TRUE, ot="Int16", dstnodata=nodata.lst[1], tile.names=tile.names, build.pyramids=FALSE)
## in parallel:
sfInit(parallel=TRUE, cpus=ifelse(length(tcovs)>24, 24, length(tcovs)))
sfExport("equi7t1", "ext", "tcovs", "mosaick.equi7t3", "make_mosaick", "tile.names", "nodata.lst", "ot.lst")
out <- sfClusterApplyLB(1:length(tcovs), function(x){try( make_mosaick(i="dominant", varn=tcovs[x], ext=ext, in.path="/data/tt/SoilGrids250m/covs1t", tr=0.002083333, r250m=TRUE, ot=ot.lst[x], dstnodata=nodata.lst[x], tile.names=tile.names, build.pyramids=FALSE) )})
sfStop()

## stack to 250m res grid
r250.lst = paste0("/data/GEOG/", c(tcovs,tcovs.DEM), "_250m_ll.tif")
detach("package:snowfall", unload=TRUE)
cl <- makeCluster(28, type="FORK")
x = parLapply(cl, 1:length(r250.lst), function(k){ system(paste0('gdalwarp ', r250.lst[k], ' /data/stacked250m/', gsub("_250m_ll.tif", ".tif", basename(r250.lst[k])), ' -co \"BIGTIFF=YES\" -r \"near\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "))) } )
stopCluster(cl)

## Factors:
system(paste0('gdalwarp /data/ESA_global/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1.tif /data/stacked250m/LCEE10.tif -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -r \"near\" -wm 2000 -overwrite -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
system(paste0('gdalwarp /data/EcoTapestry/EF_Bio_Des_250m.tif /data/stacked250m/BICUSG5.tif -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -r \"near\" -wm 2000 -overwrite -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))

