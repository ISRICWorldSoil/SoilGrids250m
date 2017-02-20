## Stacking covariates at 250 m resolution for SoilGrids
## tom.hengl@isric.org

library(RCurl)
library(rgdal)
library(raster)
library(parallel)
library(snowfall)
source("mosaick_functions.R")
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
## create new dirs:
#tile.tbl = readRDS("stacked250m_tiles.rds")
new.dirs <- paste0("/data/tt/SoilGrids250m/predicted250m/T", tile.tbl$ID)
x <- lapply(new.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

covs.lst = list.files(path="/data/stacked250m", pattern=glob2rx("*.tif$"), full.names = TRUE)
covs.lst <- covs.lst[-unlist(sapply(c("OCCGSW7","LCEE10"), function(x){grep(x, covs.lst)}))]
ov.quant <- as.list(des[match(basename(covs.lst), paste0(des$WORLDGRIDS_CODE,".tif")),"MASK_VALUE"])
names(ov.quant) = basename(covs.lst)

#sapply(c(30,32,185), function(x){ system(paste0('gdal_translate ', covs.lst[x], ' /data/tmp/', basename(covs.lst[x]), ' -ot \"Int16\" -co \"COMPRESS=DEFLATE\"'))})

## Total clean-up:
#del.lst = list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern = ".rds", recursive = TRUE, full.names = TRUE)
#unlink(del.lst)

## Function to make predictions locations
make_newdata <- function(i, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant){
  out.rds <- paste0(out.path, "/T", tile.tbl[i,"ID"], "/T", tile.tbl[i,"ID"], ".rds")
  if(!file.exists(out.rds)){
    m = readGDAL(fname=mask.l, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)
    names(m) = "LCEE10.tif"
    m$OCCGSW7.tif = readGDAL(fname=mask.w, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1
    m$OCCGSW7.tif = ifelse(m$LCEE10.tif==220|m$OCCGSW7.tif>90|m$LCEE10.tif==210, NA, m$OCCGSW7.tif)
    sel.p = !is.na(m$OCCGSW7.tif)
    if(sum(sel.p)>1){
      m = as(m, "SpatialPixelsDataFrame")
      m = m[which(sel.p),]
      x = spTransform(m, CRS("+proj=longlat +datum=WGS84"))
      m$LONWGS84 <- x@coords[,1]
      m$LATWGS84 <- x@coords[,2]
      for(j in 1:length(covs.lst)){
        m@data[,basename(covs.lst[j])] <- round(readGDAL(covs.lst[j], offset=unlist(tile.tbl[i,c("offset.y","offset.x")]), region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), silent = TRUE)$band1[m@grid.index])
      }
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
      ## Fix Landsat images (northern latitudes):
      saveRDS(m, out.rds)
      gc(); gc()
    }
  }
}
## Test:
make_newdata(i=38275, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant)
make_newdata(i=40502, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant)

## run in parallel (TAKES ca. 6 hrs):
library(snowfall)
sfInit(parallel=TRUE, cpus=54)
sfLibrary(rgdal)
sfExport("make_newdata", "tile.tbl", "ov.quant", "covs.lst")
out <- sfClusterApplyLB(as.numeric(paste(tile.tbl$ID)), function(i){ make_newdata(i, tile.tbl, covs.lst, out.path="/data/tt/SoilGrids250m/predicted250m", mask.l="/data/stacked250m/LCEE10.tif", mask.w="/data/stacked250m/OCCGSW7.tif", ov.quant) })
sfStop()
## 17532 dirs

## remove all directories with empty landmask (CAREFULL!)
pr.dirs <- basename(dirname(list.files(path="/data/tt/SoilGrids250m/predicted250m", pattern=glob2rx("*.rds$"), recursive = TRUE, full.names = TRUE)))
pr.dirs.c <- list.dirs("/data/tt/SoilGrids250m/predicted250m")[-1]
selD <- which(!basename(pr.dirs.c) %in% pr.dirs)
x = sapply(selD, function(x){unlink(pr.dirs.c[x], recursive = TRUE, force = TRUE)})
