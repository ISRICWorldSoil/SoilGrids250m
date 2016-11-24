## Fit models for USDA GreatGroups, texture classes, and several soil properties using the NASIS points (ca 350,000)
## Code by Tom.Hengl@isric.org / points prepared by Amanda Rachmaran (a.m.ramcharan@gmail.com) and Travis Nauman (tnauman@usgs.gov)

setwd("/data/NASIS")
load(".RData")
library(aqp)
library(plyr)
library(stringr)
library(dplyr)
library(sp)
library(devtools)
library(caret)
#devtools::install_github("imbs-hl/ranger/ranger-r-package/ranger") ## version to deal with Memory problems
library(ranger)
library(xgboost)
library(nnet)
library(ROCR)
library(snowfall)
library(mda)
library(psych)
library(rgdal)
library(utils)
library(R.utils)
library(raster)
library(plotKML)
library(GSIF)
library(randomForestSRC)
library(parallel)
options(rf.cores=detectCores(), mc.cores=detectCores())
source("/data/models/saveRDS_functions.R")
source("functions_NASIS.R")

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

r <- raster("/data/USA48/elev48i0100a.tif")
extent(r)
ncols = ncol(r)
nrows = nrow(r)
xllcorner = extent(r)[1]
yllcorner = extent(r)[3]
xurcorner = extent(r)[2]
yurcorner = extent(r)[4]
cellsize = res(r)[1]
NODATA_value = -32768

## legends
des <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_Covs100m.csv")
LNDCOV6.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_LandCover.csv")
PMTGSS7.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_gSSURGO_pmaterial.csv")
DRNGSS7.leg <- read.csv("/data/GEOG/NASIS/covs100m/SoilGrids_USA48_gSSURGO_drainage.csv")
PVEGKT6.leg <- read.csv("/data/GEOG/NASIS/covs100m/potential_vegetation_legend.csv")
COUNTY6.leg <- read.csv("/data/GEOG/NASIS/covs100m/Counties_legend.csv")
#GEOUSG6.leg <- read.csv("/data/GEOG/NASIS/covs100m/geology_legend.csv")
GESUSG6.leg <- read.csv("/data/GEOG/NASIS/covs100m/surfacegeology_legend.csv")

## Training points (taxonomy, texture classes):
TAX_gg.pnts = readRDS("TAX_gg.pnts.rds")
str(TAX_gg.pnts@data)
NASISpscs.pnts = readRDS("NASISpscs.pnts.rds")
NASISpscs.pnts$LOC_ID = paste("ID", NASISpscs.pnts$xwgs84, NASISpscs.pnts$ywgs84, sep="_")

## Soil properties NRCS NCSS:
system("7za e Geo_Peds07272016.csv.gz")
NCSS_peds = read.csv("Geo_Peds07272016.csv")
#summary(as.factor(NCSS_peds$geoid))
NCSS_peds$LOC_ID <- paste("ID", NCSS_peds$long, NCSS_peds$lat, sep="_")
t.props = c("clay", "sand", "soc", "bd", "n_tot") 
## clay, sand, soil organic carbon, total nitrogen, bulk density
lapply(t.props, function(i){summary(NCSS_peds[,i])})
NCSS_peds = NCSS_peds[,c("pedon_key","LOC_ID","long","lat","hzn_top","hzn_bot",t.props)]
str(NCSS_peds)

## Unique locations:
pnts = data.frame(LOC_ID=unique(c(TAX_gg.pnts$LOC_ID, NASISpscs.pnts$LOC_ID, NCSS_peds$LOC_ID)), LONWGS84=NA, LATWGS84=NA)
pnts$LONWGS84 = as.numeric(sapply(paste(pnts$LOC_ID), function(x){strsplit(x, "_")[[1]][2]}))
pnts$LATWGS84 = as.numeric(sapply(paste(pnts$LOC_ID), function(x){strsplit(x, "_")[[1]][3]}))
pnts = pnts[pnts$LONWGS84< -59 & pnts$LONWGS84> -127 & pnts$LATWGS84< 52 & pnts$LATWGS84> 24,]
str(pnts)
## 305,124 
coordinates(pnts) <- ~ LONWGS84 + LATWGS84
proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")

## Tiling system:
obj <- GDALinfo("/data/GEOG/NASIS/COUNTY6.tif")
tile.lst <- getSpatialTiles(obj, block.x=50000, return.SpatialPolygons=TRUE)
proj4string(tile.lst)
tile.tbl <- getSpatialTiles(obj, block.x=50000, return.SpatialPolygons=FALSE)
tile.tbl$ID = as.character(1:nrow(tile.tbl))
str(tile.tbl)
## 6300 tiles
tile.pol = SpatialPolygonsDataFrame(tile.lst, tile.tbl)
writeOGR(tile.pol, "tiles50km.shp", "tiles50km", "ESRI Shapefile")
plot(spTransform(pnts, CRS(proj4string(tile.lst))), pch="+")
lines(tile.lst, col="red")
## Overlay tiles and admin units (fully parallelized):
system(paste('gdal_translate /data/GEOG/NASIS/COUNTY6.tif /data/USA48/COUNTY6.sdat -ot \"Int16\" -of \"SAGA\" -a_nodata \"-32768\"'))
system(paste0(saga_cmd, ' -c=48 shapes_grid 2 -GRIDS=\"/data/USA48/COUNTY6.sgrd\" -POLYGONS=\"tiles50km.shp\" -PARALLELIZED=1 -RESULT=\"ov_COUNTY6_tiles50km.shp\"'))
ov_COUNTY6 = readOGR("ov_COUNTY6_tiles50km.shp", "ov_COUNTY6_tiles50km")
summary(sel.t <- !ov_COUNTY6@data[,"COUNTY6..MA"]==0)
## 3390 tiles with values
ov_COUNTY6 = ov_COUNTY6[sel.t,]
str(ov_COUNTY6@data)

## Prepare covariates as tiles:
t.sel = as.character(ov_COUNTY6$ID)
new.dirs <- paste0("/data/NASIS/covs100m/T", t.sel)
x <- lapply(new.dirs, dir.create, recursive=TRUE, showWarnings=FALSE)
x <- lapply(gsub("covs100m", "predicted100m", new.dirs), dir.create, recursive=TRUE, showWarnings=FALSE)
## run in parallel:
covs.lst = as.character(des$WORLDGRIDS_CODE[-which(des$WORLDGRIDS_CODE=="COUNTY6")])
#make_RDS_tiles(i=as.numeric(t.sel)[1], tile.tbl=tile.tbl, covs.lst=covs.lst, LNDCOV6.leg=LNDCOV6.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg)

#rds.lst = list.files(path="./covs100m", pattern=glob2rx("*.rds"), recursive = TRUE, full.names = TRUE)
#unlink(rds.lst)
#make_RDS_tiles(i=1922, tile.tbl=tile.tbl, covs.lst=covs.lst, LNDCOV6.leg=LNDCOV6.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg)
## TAKES 5-7 HOURS TO PREPARE
sfInit(parallel=TRUE, cpus=25)
sfExport("make_RDS_tiles", "tile.tbl", "covs.lst", "LNDCOV6.leg", "PMTGSS7.leg", "DRNGSS7.leg", "PVEGKT6.leg", "COUNTY6.leg", "GESUSG6.leg")
sfLibrary(plyr)
sfLibrary(rgdal)
sfLibrary(randomForestSRC)
sfLibrary(sp)
out <- sfClusterApplyLB(as.numeric(t.sel), function(i){ try( make_RDS_tiles(i, tile.tbl=tile.tbl, covs.lst=covs.lst, LNDCOV6.leg=LNDCOV6.leg, PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg, PVEGKT6.leg=PVEGKT6.leg, COUNTY6.leg=COUNTY6.leg, GESUSG6.leg=GESUSG6.leg) ) } )
sfStop()
save.image()

## Gap-filled Parent material map:
pm.tif = list.files(path="/data/NASIS/covs100m", pattern="PMTGSS7_predicted", recursive = TRUE, full.names = TRUE)
cat(pm.tif, sep="\n", file="my_liste.txt")
system(paste0(gdalbuildvrt, ' -input_file_list my_liste.txt PMTGSS7_f.vrt'))
system(paste0(gdalwarp, ' PMTGSS7_f.vrt PMTGSS7_f.tif -r \"near\" -co \"COMPRESS=DEFLATE\" -co \"BIGTIFF=YES\" -wm 2000 -tr ', cellsize, ' ', cellsize, ' -te ', paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")))

## OVERLAY (takes ca 15 mins):
ov <- extract.tiled(x=pnts, tile.pol=tile.pol, path="/data/NASIS/covs100m", ID="ID", cpus=48)
ov$PMTGSS7 = as.factor(ov$PMTGSS7)
ov$DRNGSS7 = as.factor(ov$DRNGSS7)
str(ov)
## "PMTGSS7" has too many levels
## http://stackoverflow.com/questions/14921805/convert-a-factor-to-indicator-variables
PMTGSS7.mat = data.frame(model.matrix(~PMTGSS7, ov))
str(PMTGSS7.mat)
PMTGSS7.mat$X.Intercept. = NULL
xsum.PMTGSS7 = sapply(PMTGSS7.mat, sum, na.rm=TRUE)
## 10 categories have <10 observations
sel.rm.PMTGSS7 = names(PMTGSS7.mat)[which(xsum.PMTGSS7<6)]
## Land cover:
summary(ov$LNDCOV6)
LNDCOV6.mat = data.frame(model.matrix(~LNDCOV6, ov))
LNDCOV6.mat$X.Intercept. = NULL
xsum.LNDCOV6 = sapply(LNDCOV6.mat, sum, na.rm=TRUE)
sel.rm.LNDCOV6 = names(LNDCOV6.mat)[which(xsum.LNDCOV6<6)]
## Potential vegetation:
PVEGKT6.mat = data.frame(model.matrix(~PVEGKT6, ov))
PVEGKT6.mat$X.Intercept. = NULL
xsum.PVEGKT6 = sapply(PVEGKT6.mat, sum, na.rm=TRUE)
sel.rm.PVEGKT6 = names(PVEGKT6.mat)[which(xsum.PVEGKT6<6)]
## Drainage class:
summary(ov$DRNGSS7)
DRNGSS7.mat = data.frame(model.matrix(~DRNGSS7, ov))
DRNGSS7.mat$X.Intercept. = NULL
xsum.DRNGSS7 = sapply(DRNGSS7.mat, sum, na.rm=TRUE)
sel.rm.DRNGSS7 = names(DRNGSS7.mat)[which(xsum.DRNGSS7<6)]

## Merge all indicators with other covariates
ov.fs = plyr::join_all(list(ov, 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(PMTGSS7.mat))]), PMTGSS7.mat[,names(PMTGSS7.mat)[which(xsum.PMTGSS7>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(LNDCOV6.mat))]), LNDCOV6.mat[,names(LNDCOV6.mat)[which(xsum.LNDCOV6>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(PVEGKT6.mat))]), PVEGKT6.mat[,names(PVEGKT6.mat)[which(xsum.PVEGKT6>5)]]), 
      cbind(data.frame(LOC_ID=ov$LOC_ID[as.numeric(rownames(DRNGSS7.mat))]), DRNGSS7.mat[,names(DRNGSS7.mat)[which(xsum.DRNGSS7>5)]])))
## join sampled values and covariates
ovA_NCSS_peds = plyr::join(NCSS_peds, ov.fs, by="LOC_ID", type="left")
ovA_TAX_gg.pnts = plyr::join(as.data.frame(TAX_gg.pnts), ov.fs, by="LOC_ID", type="left")
ovA_NASISpscs.pnts = plyr::join(as.data.frame(NASISpscs.pnts), ov.fs, by="LOC_ID", type="left")
saveRDS.gz(ovA_NCSS_peds, file="ovA_NCSS_peds.rds")
saveRDS.gz(ovA_TAX_gg.pnts, file="ovA_TAX_gg.pnts.rds")
saveRDS.gz(ovA_NASISpscs.pnts, file="ovA_NASISpscs.pnts.rds")

## covariates - WE DO NOT USE ANY FACTORS ONLY INDICATORS
pr.lst <- des$WORLDGRIDS_CODE
## clean up list of covariates
pr.lst = c(paste(pr.lst[-c(grep(pattern = "GESUSG6",pr.lst), grep(pattern = "PMTGSS7",pr.lst), grep(pattern = "PVEGKT6",pr.lst), grep(pattern = "LNDCOV6",pr.lst), grep(pattern = "DRNGSS7",pr.lst), grep(pattern = "COUNTY6",pr.lst), grep(pattern = "BLDFIE",pr.lst), grep(pattern = "ORCDRC",pr.lst), grep(pattern = "CLYPPT",pr.lst), grep(pattern = "SNDPPT",pr.lst), grep(pattern = glob2rx("N??PR5"),pr.lst), grep(pattern = glob2rx("N??PRI5"),pr.lst), grep(pattern = glob2rx("V??PRI5"),pr.lst))]), names(PMTGSS7.mat)[which(xsum.PMTGSS7>5)], names(PVEGKT6.mat)[which(xsum.PVEGKT6>5)], names(LNDCOV6.mat)[which(xsum.LNDCOV6>5)], names(DRNGSS7.mat)[which(xsum.DRNGSS7>5)])
## 205 covs

## FIT MODELS FOR FACTORS:
t.fvars = c("soiltype", "textype")
t.names = list("USDA Great Group", "SPCS classes")
names(t.names) = t.fvars
ov.lst = list(ovA_TAX_gg.pnts, ovA_NASISpscs.pnts)
names(ov.lst) = t.fvars
save.image()
## fit models in a loop:
for(j in t.fvars){
  if(!file.exists(paste0("mRF_", j,".rds"))){
    formulaString = as.formula(paste(j, ' ~ ', paste(pr.lst, collapse="+")))
    x = ov.lst[[j]][,all.vars(formulaString)] ## sample.int(nrow(ov.lst[[j]]), 2e4)
    x[,j] = droplevels(x[,j])
    x = x[complete.cases(x),]
    mRF = rfsrc(formulaString, data=x, ntree = 500, membership=TRUE) ## ov.lst[[j]][sample.int(nrow(ov.lst[[j]]), 2e4),all.vars(formulaString)]
    #mRF = ranger::ranger(formulaString, data=x, importance="impurity", write.forest=TRUE, probability=TRUE) ## respect.unordered.factors=TRUE
    #save(formulaString, x, file="classification.RData")
    ## https://www.r-bloggers.com/on-ranger-respect-unordered-factors/
    gc(); gc()
    cat("Results of model fitting 'randomForest':\n", file=paste0("NASIS_resultsFit_ ", j, ".txt"))
    cat("\n", file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    cat(paste("Variable:", t.names[[j]]), file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    cat("\n", file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    sink(file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE, type="output")
    cat("\n Random forest model:", file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    print(mRF)
    cat("\n Variable importance:\n", file=paste0("NASIS_resultsFit_ ", j, ".txt"), append=TRUE)
    s.mRF = rfsrc(formulaString, data=x[sample.int(nrow(x), 5e4),], ntree = 500)
    xl <- vimp(s.mRF)
    print(data.frame(xl$importance)[,1], max.levels=length(all.vars(formulaString)))
    write.csv(data.frame(xl$importance), file=paste0("vimp_mRF_", j,".csv"))
    #xl <- as.list(ranger::importance(mRF))
    #print(t(data.frame(xl[order(unlist(xl), decreasing=TRUE)[1:25]])))
    sink()
    saveRDS.gz(mRF, file=paste0("mRF_", j,".rds"))
  }
}

## Predict classes in parallel:
mRF = readRDS.gz("mRF_soiltype.rds")
#str(mRF)
levs = attr(mRF$predicted, "dimnames")[[2]]
varn = "TAXgg"
#predict_factor_tile(i=1042, mRF, varn="TAXgg", levs=levs, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg)
gc(); gc()
save.image()

sfInit(parallel=TRUE, cpus=6)
sfExport("predict_factor_tile", "mRF", "PMTGSS7.leg", "DRNGSS7.leg", "levs", "t.sel", "varn")
sfLibrary(plyr)
sfLibrary(rgdal)
sfLibrary(randomForestSRC)
out <- sfClusterApplyLB(as.numeric(t.sel), function(i){ try( predict_factor_tile(i, mRF, varn=varn, levs=levs, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg) ) } )
sfStop()

## Make mosaics:
te = paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")
sfInit(parallel=TRUE, cpus=ifelse(length(levs)>46, 46, length(levs)))
sfExport("gdalbuildvrt", "gdalwarp", "levs", "mosaic_tiles_100m", "varn", "te")
out <- sfClusterApplyLB(levs, function(x){try( mosaic_tiles_100m(x, in.path="/data/NASIS/predicted100m", varn=varn, out.path="/data/GEOG/NASIS/predicted100m/", te=te) )})
sfStop()
## most probable class:
tmp.lst <- c(list.files(path="/data/NASIS/predicted100m", pattern=glob2rx(paste0(varn, "_T????.tif$")), full.names=TRUE, recursive=TRUE), list.files(path="/data/NASIS/predicted100m", pattern=glob2rx(paste0(varn, "_T???.tif$")), full.names=TRUE, recursive=TRUE))
out.tmp <- tempfile(fileext = ".txt")
vrt.tmp <- tempfile(fileext = ".vrt")
cat(tmp.lst, sep="\n", file=out.tmp)
system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
system(paste0(gdalwarp, ' ', vrt.tmp, ' /data/GEOG/NASIS/predicted100m/TAXgg_100m.tif -ot \"Int16\" -dstnodata \"32767\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -tr 100 100 -r \"near\" -te ', te))
write.csv(data.frame(Value=1:length(levs), Class=levs), file="/data/GEOG/NASIS/predicted100m/TAXgg_100m.tif.csv")
save.image()

## Soil texture families:
mRF = readRDS.gz("mRF_textype.rds")
#str(mRF)
levs = attr(mRF$predicted, "dimnames")[[2]]
varn = "PSCS"
#predict_factor_tile(i=1042, mRF, varn=varn, levs=levs, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg)
gc(); gc()

## Because textype has only 20+ classes we can use more cores
sfInit(parallel=TRUE, cpus=12)
sfExport("predict_factor_tile", "mRF", "PMTGSS7.leg", "DRNGSS7.leg", "levs", "t.sel", "varn")
sfLibrary(plyr)
sfLibrary(rgdal)
sfLibrary(randomForestSRC)
out <- sfClusterApplyLB(as.numeric(t.sel), function(i){ try( predict_factor_tile(i, mRF, varn=varn, levs=levs, in.path="/data/NASIS/covs100m", out.path="/data/NASIS/predicted100m", PMTGSS7.leg=PMTGSS7.leg, DRNGSS7.leg=DRNGSS7.leg) ) } )
sfStop()
save.image()

## Make mosaics:
te = paste(as.vector(extent(r))[c(1,3,2,4)], collapse=" ")
levsF = basename(list.files(path="/data/NASIS/predicted100m/T1042", pattern=glob2rx("PSCS_*_T1042.tif$")))
levsF = sapply(levsF, function(x){strsplit(x, "_")[[1]][2]})
sfInit(parallel=TRUE, cpus=ifelse(length(levsF)>46, 46, length(levsF)))
sfExport("gdalbuildvrt", "gdalwarp", "levsF", "mosaic_tiles_100m", "varn", "te")
out <- sfClusterApplyLB(levsF, function(x){try( mosaic_tiles_100m(x, in.path="/data/NASIS/predicted100m", varn=varn, out.path="/data/GEOG/NASIS/predicted100m/", te=te) )})
sfStop()
## most probable class:
tmp.lst <- c(list.files(path="/data/NASIS/predicted100m", pattern=glob2rx(paste0(varn, "_T????.tif$")), full.names=TRUE, recursive=TRUE), list.files(path="/data/NASIS/predicted100m", pattern=glob2rx(paste0(varn, "_T???.tif$")), full.names=TRUE, recursive=TRUE))
out.tmp <- tempfile(fileext = ".txt")
vrt.tmp <- tempfile(fileext = ".vrt")
cat(tmp.lst, sep="\n", file=out.tmp)
system(paste0(gdalbuildvrt, ' -input_file_list ', out.tmp, ' ', vrt.tmp))
system(paste0(gdalwarp, ' ', vrt.tmp, ' /data/GEOG/NASIS/predicted100m/PSCS_100m.tif -ot \"Int16\" -dstnodata \"32767\" -co \"BIGTIFF=YES\" -multi -wm 2000 -co \"COMPRESS=DEFLATE\" -tr 100 100 -r \"near\" -te ', te))
write.csv(data.frame(Value=1:length(levs), Class=levs), file="/data/GEOG/NASIS/predicted100m/PSCS_100m.tif.csv")
save.image()


## FIT MODELS FOR NUMERIC VARIABLES:



