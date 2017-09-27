## Stacking covariates at 250 m resolution for SoilGrids
## tom.hengl@isric.org

list.of.packages <- c("RCurl", "rgdal", "raster", "parallel", "snowfall", "plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

setwd("/data/models")
load(".RData")
library(RCurl)
library(rgdal)
library(raster)
library(parallel)
library(snowfall)
source("mosaick_functions.R")
## Grid definition:
r <- raster("/data/stacked250m/LCEE10.tif")
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
des.s = des[!des$ROOT_FILE=="",c("WORLDGRIDS_CODE","ROOT_FILE","RESAMPLE_METHOD","DATA_FORMAT")]
s = file.exists(paste(des.s$ROOT_FILE))
summary(s)
#des.s$ROOT_FILE[!s]
#system(paste0('gdalinfo /data/MOD13Q1/M_listeJanFeb.vrt'))
save.image()

## process in parallel (250m):
#ERROR 2: CPLMalloc(): Out of memory allocating 90720000 bytes.
sfInit(parallel=TRUE, cpus=20)
sfExport("r", "cellsize", "des.s")
sfLibrary(rgdal)
sfLibrary(raster)
out <- sfClusterApplyLB(1:nrow(des.s), function(k){ if(!file.exists(paste0('/data/stacked250m/', des.s$WORLDGRIDS_CODE[k], '.tif'))){ system(paste0('gdalwarp ', des.s$ROOT_FILE[k], ' /data/stacked250m/', des.s$WORLDGRIDS_CODE[k], '.tif -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -r \"', des.s$RESAMPLE_METHOD[k], '\" -ot \"', des.s$DATA_FORMAT[k], '\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "))) } } )
sfStop()

## prepare covariates in 1 km resolution
# sfInit(parallel=TRUE, cpus=ifelse(length(tcovs)>24, 24, length(tcovs)))
# sfExport("r", "des.s")
# sfLibrary(rgdal)
# sfLibrary(raster)
# out <- sfClusterApplyLB(1:nrow(des.s), function(k){ if(!file.exists(paste0('/data/stacked1km/', des.s$WORLDGRIDS_CODE[k], '.tif'))){ system(paste0('gdalwarp ', des.s$ROOT_FILE[k], ' /data/stacked1km/', des.s$WORLDGRIDS_CODE[k], '.tif -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', 1/120, ' ', 1/120, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "))) } } )
# sfStop()

## Covariates in Equi7 system (tiles)
# library(snowfall)
# sel.equi7 = c(grep(glob2rx("***USG5"), des$WORLDGRIDS_CODE), grep(glob2rx("***MRG5"), des$WORLDGRIDS_CODE))
# tcovs <- as.character(des$WORLDGRIDS_CODE[sel.equi7])
# nodata.lst <- des$NO_DATA[sel.equi7]
# ot.lst <- as.character(des$DATA_FORMAT[sel.equi7])
## test it:
#make_mosaick(i="dominant", varn=tcovs[1], ext=ext, in.path="/data/tt/SoilGrids250m/covs1t", tr=0.002083333, r250m=TRUE, ot="Int16", dstnodata=nodata.lst[1], tile.names=tile.names, build.pyramids=FALSE)
## in parallel:
# sfInit(parallel=TRUE, cpus=ifelse(length(tcovs)>24, 24, length(tcovs)))
# sfExport("equi7t1", "ext", "tcovs", "mosaick.equi7t3", "make_mosaick", "tile.names", "nodata.lst", "ot.lst")
# out <- sfClusterApplyLB(1:length(tcovs), function(x){try( make_mosaick(i="dominant", varn=tcovs[x], ext=ext, in.path="/data/tt/SoilGrids250m/covs1t", tr=0.002083333, r250m=TRUE, ot=ot.lst[x], dstnodata=nodata.lst[x], tile.names=tile.names, build.pyramids=FALSE) )})
# sfStop()

## stack to 250m res grid
r250.lst = list.files(path="/data/GEOG", pattern=glob2rx("???MRG5_250m_ll.tif"), full.names=TRUE)
detach("package:snowfall", unload=TRUE)
cl <- parallel::makeCluster(length(r250.lst), type="FORK")
x = parallel::parLapply(cl, 1:length(r250.lst), function(k){ system(paste0('gdalwarp ', r250.lst[k], ' /data/stacked250m/', gsub("_250m_ll.tif", ".tif", basename(r250.lst[k])), ' -co \"BIGTIFF=YES\" -r \"near\" -wm 2000 -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" "))) } )
stopCluster(cl)

## Factors:
system(paste0('gdalwarp /data/ESA_global/ESACCI-LC-L4-LCCS-Map-300m-P5Y-2010-v1.6.1.tif /data/stacked250m/LCEE10.tif -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -r \"near\" -wm 2000 -overwrite -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))
system(paste0('gdalwarp /data/EcoTapestry/EF_Bio_Des_250m.tif /data/stacked250m/BICUSG5.tif -t_srs \"', proj4string(r), '\" -co \"BIGTIFF=YES\" -r \"near\" -wm 2000 -overwrite -co \"COMPRESS=DEFLATE\" -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))

## Tiling system:
obj <- GDALinfo("/data/stacked250m/LCEE10.tif")
tile.lst <- getSpatialTiles(obj, block.x=1, return.SpatialPolygons=TRUE)
tile.tbl <- getSpatialTiles(obj, block.x=1, return.SpatialPolygons=FALSE)
tile.tbl$ID = as.character(1:nrow(tile.tbl))
str(tile.tbl)
## 54000 tiles
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
unlink("tiles_ll_100km.shp")
writeOGR(tile.pol, "tiles_ll_100km.shp", "tiles_ll_100km", "ESRI Shapefile")
saveRDS(tile.tbl, "stacked250m_tiles.rds")
saveRDS(tile.pol, "stacked250m_tiles_pol.rds")
## create new dirs:
#tile.tbl = readRDS("stacked250m_tiles.rds")
new.dirs <- paste0("/data/tt/SoilGrids250m/predicted250m/T", tile.tbl$ID)
x <- lapply(new.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

covs.lst = list.files(path="/data/stacked250m", pattern=glob2rx("*.tif$"), full.names = TRUE)
## 218
## Check if all layers stack:
#ts = stack(covs.lst)
#ext.lst = lapply(covs.lst, function(i){paste(as.vector(raster::extent(raster(i))), collapse="_")})
#del.lst = covs.lst[sapply(ext.lst, function(i){!i==paste(as.vector(raster::extent(r)), collapse="_")})]
## subsample values:
#s.pnt = sampleRandom(r, size=200, sp=TRUE, na.rm=TRUE)
#xs = mclapply(covs.lst, FUN=function(i){raster::extract(raster(i), s.pnt)}, mc.cores = 46)
#names(xs) = covs.lst
#xs = as.data.frame(do.call(rbind, lapply(xs, summary)))
#saveRDS(xs, file="covs_statistics.rds")
xs = readRDS("covs_statistics.rds")
#write.csv(xs, "covs_statistics.csv")

summary((basename(covs.lst)) %in% paste0(des$WORLDGRIDS_CODE, ".tif"))
## clean-up:
unlink(covs.lst[!((basename(covs.lst)) %in% paste0(des$WORLDGRIDS_CODE, ".tif"))])
covs.lst = list.files(path="/data/stacked250m", pattern=glob2rx("*.tif$"), full.names = TRUE)
covs.lst <- covs.lst[-unlist(sapply(c("N11MSD3","LCEE10","CSCMCF5","B02CHE3","B08CHE3","B09CHE3","S01ESA4","S02ESA4","S11ESA4","S12ESA4"), function(x){grep(x, covs.lst)}))]
str(covs.lst)
## 208
#ov.quant <- as.list(des[match(basename(covs.lst), paste0(des$WORLDGRIDS_CODE,".tif")),"MASK_VALUE"])
ov.quant <- as.list(xs$Mean[which(list.files(path="/data/stacked250m", pattern=glob2rx("*.tif$"), full.names = TRUE) %in% covs.lst)])
names(ov.quant) = basename(covs.lst)

#sapply(c(30,32,185), function(x){ system(paste0('gdal_translate ', covs.lst[x], ' /data/tmp/', basename(covs.lst[x]), ' -ot \"Int16\" -co \"COMPRESS=DEFLATE\"'))})

## Total clean-up:
#del.lst = list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern = ".rds", recursive = TRUE, full.names = TRUE)
#unlink(del.lst)
#del.lst = list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern = ".tif", recursive = TRUE, full.names = TRUE)
#unlink(del.lst)

## Function to make predictions locations ----
make_newdata <- function(i, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant){
  out.rds <- paste0(out.path, "/T", tile.tbl[i,"ID"], "/T", tile.tbl[i,"ID"], ".rds")
  if(!file.exists(out.rds)){
    m = readGDAL(fname=mask.l, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    names(m) = "LCEE10.tif"
    m$OCCGSW7.tif = readGDAL(fname=mask.w, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1
    ## DEFINITION OF LAND MASK:
    #m$OCCGSW7.tif = ifelse(m$LCEE10.tif==220|m$OCCGSW7.tif>95|m$LCEE10.tif==210, NA, m$OCCGSW7.tif)
    m$LCEE10.tif = ifelse(m$LCEE10.tif==220|m$LCEE10.tif==210, NA, m$LCEE10.tif)
    sel.p = !is.na(m$LCEE10.tif)
    if(sum(sel.p)>1){
      m = as(m, "SpatialPixelsDataFrame")
      m = m[which(sel.p),]
      x = spTransform(m, CRS("+proj=longlat +datum=WGS84"))
      m$LONWGS84 <- x@coords[,1]
      m$LATWGS84 <- x@coords[,2]
      for(j in 1:length(covs.lst)){
        m@data[,basename(covs.lst[j])] <- round(readGDAL(covs.lst[j], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[m@grid.index])
      }
      ## Surface water images (https://global-surface-water.appspot.com/download):
      ## OCCGSW7.tif == 101 => missing values
      m$CHAGSW7.tif = ifelse(m$CHAGSW7.tif>200, 100, m$CHAGSW7.tif)
      m$EXTGSW7.tif = ifelse(m$EXTGSW7.tif==1, 1, 0)
      m$OCCGSW7.tif = ifelse(m$OCCGSW7.tif>100, 0, m$OCCGSW7.tif)
      ## Fill-in missing values (if necessary):
      sel.mis = sapply(m@data, function(x){sum(is.na(x))>0})
      if(sum(sel.mis)>0){
        for(j in which(sel.mis)){
          if(!is.factor(m@data[,j])){
            repn = quantile(m@data[,j], probs=.5, na.rm=TRUE)
            if(is.na(repn)){
              repn = ov.quant[[names(m)[j]]]
            }
            m@data[,j] = ifelse(is.na(m@data[,j]), repn, m@data[,j])
          }
        }
      }
      ## Fix Landsat images (northern latitudes)?
      ## REDL00, NIRL00, SW1L00, SW2L00 == 0 => missing value
      saveRDS(m, out.rds)
      gc(); gc()
    }
  } else {
    ## check if maybe some layer is missing:
    m <- readRDS(out.rds)
    mis.cov = covs.lst[which(!basename(covs.lst) %in% names(m))]
    if(length(mis.cov)>0){
      for(k in 1:length(mis.cov)){
        m@data[,basename(mis.cov[k])] <- round(readGDAL(mis.cov[k], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[m@grid.index])
      }
      saveRDS(m, out.rds)
    }
  }
}
save.image()

## Test:
#make_newdata(i=38275, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant)
#make_newdata(i=40502, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant)
#make_newdata(i=45489, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant)
## Africa border with Sahara:
#make_newdata(i=27198, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant)

## run in parallel (TAKES ca. 28 hrs!!):
library(snowfall)
snowfall::sfInit(parallel=TRUE, cpus=48)
snowfall::sfLibrary(rgdal)
snowfall::sfExport("make_newdata", "tile.tbl", "ov.quant", "covs.lst")
out <- snowfall::sfClusterApplyLB(as.numeric(paste(tile.tbl$ID)), function(i){ make_newdata(i, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant) })
#out <- snowfall::sfClusterApplyLB(as.numeric(sapply(pr.dirs, function(j){strsplit(j, "T")[[1]][2]})), function(i){ make_newdata(i, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant) })
snowfall::sfStop()

## remove all directories with empty landmask (CAREFULL - only run once all tiles are finished!)
pr.dirs <- basename(dirname(list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
pr.dirs.c <- list.dirs("/data/tt/SoilGrids250m/predicted250m")[-1]
selD <- which(!basename(pr.dirs.c) %in% pr.dirs)
x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})
saveRDS(pr.dirs, "prediction_dirs.rds")
## 18,653 dirs on the end